#include "placer.h"

void Placer::run()
{
    derived = DERIVED_PLACE;

    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "Added reads\t#Paths:" << added_reads->size() << "\n";
        PathReadsMap::iterator jt;
        size_t i;
        for ( i = 1, jt = added_reads->begin(); jt != added_reads->end(); ++i, ++jt ) {
            std::cout << i << "\tPathId:" << jt->first << "\t";
            std::cout << "# reads:" << jt->second.size() << "\n";
            ReadEntryList::iterator rt;
            for ( rt = jt->second.begin(); rt != jt->second.end(); ++rt )
                std::cout << rt->read << "\t" << rt->rpos << "\t" << rt->ppos << "\n";
        }
    }

    ct_gappy = ct_gappy_merge = ct_gappy_extend = 0;

    ct_realign = ct_empty_ref = 0;
    ct_msa_empty = ct_msa_lerr = ct_msa_rerr = ct_msa_hole = ct_msa_ltrim = ct_msa_rtrim = 0;
    ct_con_lerr = ct_con_rerr = ct_con_hole = ct_con_ltrim = ct_con_rtrim = 0;

    align();
    check();
    write();

    status = true;
}

void Placer::align()
{
    progress = Progress( 1.0, 0, npaths, mytime() );

    size_t i;
#pragma omp parallel for schedule(dynamic, 1) if (Param::ncpus>1) private(i) num_threads(Param::ncpus)
    for ( i = 0; i < npaths; i++ ) {
        assert( Id2Paths.find(i) != Id2Paths.end() );
        PathId pid = Id2Paths[i];
        size_t nid = Path2Ids[pid];
        //assert(i == nid);
        alignAugmentPath( nid, Id2Paths[i], Alns[i] );

        double prev = progress.ratio;
#pragma omp critical
        {
            progress.count++;
            progress.showProgress();
            if ( progress.ratio > prev ) {
                if (Param::verbose >= 1 && ct_realign > 0 ) 
                    printf("# re-aligned paths:%zu\n", ct_realign);
                if (Param::verbose > 1 && ct_gappy > 0 ) 
                    printf( "Gapped placements:%zu (clustering:%zu, extension:%zu)\n", ct_gappy, ct_gappy_merge, ct_gappy_extend );
            }
        }
    }

    printf("# re-aligned paths:%zu\n", ct_realign);

    int assem = countAssembledReads();
    double ratio = (double)assem/nreads;
    printf("assembled read ratio: assembled/total-reads = %.2f%% (%d/%d)\n", ratio*100, assem, nreads);
}


AlignSummary Placer::findAlignSummary( PathId pid, PathAlignMap &aligns)
{
    AlignSummary sum;
    PathAlignMap::iterator it = aligns.find(pid);

    if ( it != aligns.end() ) {
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "AlignSummary found for pid:" << pid << "\n";
        sum = it->second;
    } else {
        if ( Param::verbose > 1 || Param::debug_flag )
            std::cout << "AlignSummary NOT found for pid:" << pid << "\n";
    }
    return sum;
}


void Placer::alignAugmentPath( size_t nid, PathId lid, PathAligner &paln )
{
    if ( Param::debug_id != -1 && Param::debug_id != lid ) return;

    // Param::extend_first ? 
    //     alignFinalMergers(nid,lid,paln) :
    //     alignFinalJoiners(nid,lid,paln) ;
    if ( Param::extend_first ) 
        alignFinalMergers(nid,lid,paln);
    else {
        Param::read_bridge_extend ? 
            alignReadsJoiners(nid,lid,paln) :
            alignFinalJoiners(nid,lid,paln) ;
    }
}

void Placer::alignReadsJoiners( size_t nid, PathId lid, PathAligner &paln )
{
    if ( Param::debug_id != -1 && Param::debug_id != lid ) return;

    ExtenderMap::iterator mt = read_exts->find(lid);
    assert( mt != read_exts->end() );
    Extender ex = mt->second;

    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "\nExtended path:" << "s" << lid << "\n";
        std::cout << ex.sequence << "\n";
        std::cout << "#Members:" << ex.members.size() << "\t";
        for ( PathIdSet::iterator pt = ex.members.begin(); 
              pt != ex.members.end(); ++pt ) 
            std::cout << "s" << *pt << ":" << ex.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << ex.aligns.size() << "\n";
    }
        
    paln.setSequence(ex.sequence);
    paln.setLStop(ex.lstop);
    paln.setRStop(ex.rstop);
    paln.setLTrim(ex.ltrim);
    paln.setRTrim(ex.rtrim);

    bool gappy = false;
    std::multimap<int, PathId> posmap = util::sortByValue<PathId,int>(ex.positions);
    std::multimap<int, PathId>::iterator it;

    int k = 0;
    for ( it = posmap.begin(); it != posmap.end(); ++it ) { 
        PathId pid = it->second;
        int    pos = it->first;
        AlignSummary sum = findAlignSummary(pid, ex.aligns);
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Member:s" << pid << "\tpos:" << pos << "\n";

        assert( read_idmap->find(pid) != read_idmap->end() );
        PathId epid = (*read_idmap)[pid];
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Short-pid:s" << epid << "\n";

        PathAligner ealn;
        alignShortExtensionCluster(epid, gappy, ealn);

        ealn.adjust(pos);
        paln.merge(ealn);

        k++;
    }

    if ( added_reads->find(lid) != added_reads->end() ) {
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Additional reads:" << (*added_reads)[lid].size() << "\n";
        ReadEntryList::iterator rt;
        size_t i = 0;
        for ( rt = (*added_reads)[lid].begin(); rt != (*added_reads)[lid].end(); ++rt ) {
            ReadPlacement place;
            place.rid = rt->read;
            place.read_pos = rt->rpos;
            place.path_pos = rt->ppos;
            
            if ( Param::verbose > 1 ) {
                std::cout << ++i << "\t";
                place.print(std::cout);
            }
            paln.addPlacement(place);
        }
    }

    if ( gappy ){
#pragma omp atomic
        ct_gappy++;
        
        //if ( Param::summary_flag ) {
        if ( Param::verbose > 1 || Param::debug_flag ) {
#pragma omp critical
            std::cout << "Id:" << nid << "\tPath:s" << lid << "\tgappy alignment\n";
        }
        ct_realign++;
        paln.realign(seqs);
    }

    if ( Param::verbose > 1 ) {
        std::cout << "Extend path placement\n";
        paln.printPlacement(std::cout);
    }
}
    
/*
  ===========================oo J1
  xxxxxxxxxxxxooooo C1
              Long overlap
              oooooxxxxxxxxxxx C2
                           Short overlap beteween J1 & J2
                             oo================================= J2
                             xxxxxxxxxxxxxxoooo C3
                                           Long overalap bwteen C3 & C4
                                           oooooxxxxxxxxxxx C4 (Cluster: P1,P2,P3)
                                           ---------------- P1
                                             ------         P2
                                                 ---------  P3
 */

void Placer::alignFinalJoiners( size_t nid, PathId lid, PathAligner &paln )
{
    if ( Param::debug_id != -1 && Param::debug_id != lid ) return;

    ExtenderMap::iterator mt = tiny_exts->find(lid);
    assert( mt != tiny_exts->end() );
    Extender ex = mt->second;

    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "\nExtended path:" << "s" << lid << "\n";
        std::cout << ex.sequence << "\n";
        std::cout << "#Members:" << ex.members.size() << "\t";
        for ( PathIdSet::iterator pt = ex.members.begin(); 
              pt != ex.members.end(); ++pt ) 
            std::cout << "s" << *pt << ":" << ex.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << ex.aligns.size() << "\n";
    }
        
    paln.setSequence(ex.sequence);
    paln.setLStop(ex.lstop);
    paln.setRStop(ex.rstop);
    paln.setLTrim(ex.ltrim);
    paln.setRTrim(ex.rtrim);

    bool gappy = false;
    std::multimap<int, PathId> posmap = util::sortByValue<PathId,int>(ex.positions);
    std::multimap<int, PathId>::iterator it;

    int k = 0;
    for ( it = posmap.begin(); it != posmap.end(); ++it ) { 
        PathId pid = it->second;
        int    pos = it->first;
        AlignSummary sum = findAlignSummary(pid, ex.aligns);
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Member:s" << pid << "\tpos:" << pos << "\n";

        assert( tiny_idmap->find(pid) != tiny_idmap->end() );
        PathId epid = (*tiny_idmap)[pid];
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Long-pid:l" << epid << "\n";

        PathAligner ealn;
        alignLongExtensionCluster(epid, gappy, ealn);

        if ( sum.dlist.size() || sum.ilist.size() ) {
            ct_gappy_extend++;
            gappy = true;
        }

        ealn.adjust(pos);
        paln.merge(ealn);

        k++;
    }

    if ( gappy ){
#pragma omp atomic
        ct_gappy++;
        
        //if ( Param::summary_flag ) {
        if ( Param::verbose > 1 || Param::debug_flag ) {
#pragma omp critical
            std::cout << "Id:" << nid << "\tPath:s" << lid << "\tgappy alignment\n";
        }
        ct_realign++;
        paln.realign(seqs);
    }

    if ( Param::verbose > 1 ) {
        std::cout << "Extend path placement\n";
        paln.printPlacement(std::cout);
    }

}

void Placer::alignFinalMergers( size_t nid, PathId cpivot, PathAligner &paln )
{
    if ( Param::debug_id != -1 && Param::debug_id != cpivot ) return;

    alignSuperCluster(nid, cpivot, paln);
}

void Placer::alignSuperCluster( size_t nid, PathId cpivot, PathAligner &paln )
{
    Clusters::iterator ct = clusters->find(cpivot);
    assert( ct != clusters->end() );

    Cluster cluster = ct->second;
    //assert( path_entries->find(cpivot) != path_entries->end() );
    std::string cpivot_seq = cluster.sequence;
    assert(cpivot_seq.size()>0);

    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "\nSuper cluster:c" << cpivot << "\n";
        std::cout << cpivot_seq << "\n";
        std::cout << "#Members:" << cluster.members.size() << "\t";
        PathIdList::iterator pt;
        AlignSummaryList::iterator at;
        for ( pt  = cluster.members.begin(), at  = cluster.summarys.begin(); 
              pt != cluster.members.end()  , at != cluster.summarys.end(); 
              ++pt, ++at ) {
            std::cout << "c" << *pt << ":" << at->lgap.second << " ";
        }
        std::cout << "\n";
    }

    assert( cluster_idmap->find(cpivot) != cluster_idmap->end() );
    //PathId opid = (*cluster_idmap)[cpivot];
    paln.setSequence(cpivot_seq);
    paln.setLStop( cluster.lstop );
    paln.setRStop( cluster.rstop );
    paln.setLTrim( cluster.ltrim );
    paln.setRTrim( cluster.rtrim );

    bool gappy = false;
    PathIdList::iterator pt;
    AlignSummaryList::iterator at;
    for ( pt  = cluster.members.begin(), at  = cluster.summarys.begin(); 
          pt != cluster.members.end()  , at != cluster.summarys.end(); 
          ++pt, ++at ) {
        PathId cmember = *pt;
        AlignSummary csum = *at;
        int cpos = csum.lgap.second;
        assert( cpos >= 0 );

        if ( Param::verbose > 1 ) csum.print(std::cout);
        if ( Param::verbose > 1 || Param::debug_flag ) {
            std::cout << "Cluster:c" << cmember << "\tpos:" << cpos << "\n";
            std::cout << cluster.sequence << "\n";
        }

        if ( csum.ilist.size() || csum.dlist.size() ) {
            ct_gappy_merge++;
            gappy = true;
        }
        if ( Param::verbose > 1 && gappy ) std::cout << "gappy\n";

        assert ( cluster_idmap->find(cmember) != cluster_idmap->end() );
        PathId tpivot = (*cluster_idmap)[cmember];
        PathAligner taln;
        if ( !Param::read_bridge_extend ) 
            alignTinyExtensionPath( tpivot, taln, gappy );
        else
            alignReadExtensionPath( tpivot, taln, gappy );

        taln.adjust(cpos);
        paln.merge(taln);
        
    }

    if ( gappy ){
#pragma omp atomic
        ct_gappy++;

        if ( Param::verbose > 1 || Param::debug_flag ) {
#pragma omp critical
            std::cout << "Id:" << nid << "\tPath:c" << cpivot << "\tgappy alignment\n";
        }
        ct_realign++;
        paln.realign(seqs);
    }

    if ( Param::verbose > 1 ) {
        std::cout << "Super cluster placement\n";
        paln.printPlacement(std::cout);
    }
}

void Placer::alignReadExtensionPath( PathId tpivot,
                                     PathAligner &taln,
                                     bool &gappy )
{
    setbuf(stdout, NULL);

    if ( Param::verbose > 1 || Param::debug_flag ) 
        std::cout << "Bridge-pivot:r" << tpivot << "\n";
    ExtenderMap::iterator tit = read_exts->find(tpivot);
    assert( tit != read_exts->end() );
    
    Extender ext = tit->second;
    std::multimap<int, PathId> tposs = util::sortByValue<PathId,int>(ext.positions);
    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "#Members:" << ext.members.size() << "\n";
        for ( PathIdSet::iterator pt = ext.members.begin(); 
              pt != ext.members.end(); ++pt ) 
            std::cout << "r" << *pt << ":" << ext.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << ext.aligns.size() << "\n";
    }
    
    std::multimap<int, PathId>::iterator pit;
    for ( pit = tposs.begin(); pit != tposs.end(); ++pit ) {
        PathId tmember = pit->second;
        int    tpos = pit->first;
        assert( tpos >= 0 );
        
        if ( Param::verbose > 1 || Param::debug_flag ) {
            std::cout << "Bridge-pid:s" << tmember << "\tpos:" << tpos << "\n";
            std::cout << (*read_exts)[tmember].sequence << "\n";
        }
        AlignSummary tsum = findAlignSummary(tmember, ext.aligns);
        if ( Param::verbose > 1 ) tsum.print(std::cout);
        if ( tsum.ilist.size() || tsum.dlist.size() ) {
            ct_gappy_extend++;
            gappy = true;
        }
        if ( Param::verbose > 1 && gappy ) std::cout << "gappy\n";

        assert( read_idmap->find(tmember) != read_idmap->end() );
        PathId spivot = (*read_idmap)[tmember];
        
        PathAligner saln;
        alignTinyExtensionPath( spivot, saln, gappy );
        saln.adjust( tpos );
        taln.merge( saln );
    }

    if ( added_reads->find(tpivot) != added_reads->end() ) {
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "PathId:" << tpivot << "\tAdditional reads:" << (*added_reads)[tpivot].size() << "\n";
        ReadEntryList::iterator rt;
        size_t i = 0;
        for ( rt = (*added_reads)[tpivot].begin(); rt != (*added_reads)[tpivot].end(); ++rt ) {
            ReadPlacement place;
            place.rid = rt->read;
            place.read_pos = rt->rpos;
            place.path_pos = rt->ppos;
            
            if ( Param::verbose > 1 ) {
                std::cout << ++i << "\t";
                place.print(std::cout);
            }
            taln.addPlacement(place);
        }
    }

    if ( Param::verbose > 1 ) {
        std::cout << "Read bridge path extension placement\n";
        taln.printPlacement(std::cout);        
    }
}

void Placer::alignTinyExtensionPath( PathId tpivot,
                                     PathAligner &taln,
                                     bool &gappy )
{
    if ( Param::verbose > 1 || Param::debug_flag ) 
        std::cout << "Tiny-pivot:s" << tpivot << "\n";
    ExtenderMap::iterator tit = tiny_exts->find(tpivot);
    assert( tit != tiny_exts->end() );
    
    Extender ext = tit->second;
    std::multimap<int, PathId> tposs = util::sortByValue<PathId,int>(ext.positions);
    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "#Members:" << ext.members.size() << "\n";
        for ( PathIdSet::iterator pt = ext.members.begin(); 
              pt != ext.members.end(); ++pt ) 
            std::cout << "s" << *pt << ":" << ext.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << ext.aligns.size() << "\n";
    }
    
    std::multimap<int, PathId>::iterator pit;
    for ( pit = tposs.begin(); pit != tposs.end(); ++pit ) {
        PathId tmember = pit->second;
        int    tpos = pit->first;
        assert( tpos >= 0 );
        
        if ( Param::verbose > 1 || Param::debug_flag ) {
            std::cout << "Short-pid:s" << tmember << "\tpos:" << tpos << "\n";
            std::cout << (*tiny_exts)[tmember].sequence << "\n";
        }
        AlignSummary tsum = findAlignSummary(tmember, ext.aligns);
        if ( Param::verbose > 1 ) tsum.print(std::cout);
        if ( tsum.ilist.size() || tsum.dlist.size() ) {
            ct_gappy_extend++;
            gappy = true;
        }
        if ( Param::verbose > 1 && gappy ) std::cout << "gappy\n";

        assert( tiny_idmap->find(tmember) != tiny_idmap->end() );
        PathId lpivot = (*tiny_idmap)[tmember];
        
        PathAligner laln;
        alignLongExtensionPath( lpivot, laln, gappy );
        laln.adjust( tpos );
        taln.merge( laln );
    }
    if ( Param::verbose > 1 ) {
        std::cout << "Short extension placement\n";
        taln.printPlacement(std::cout);        
    }
}

void Placer::alignLongExtensionPath( PathId lpivot,
                                     PathAligner &laln,
                                     bool &gappy ) 
{
    if ( Param::verbose > 1 || Param::debug_flag ) 
        std::cout << "Long-pivot:l" << lpivot << "\n";

    ExtenderMap::iterator lxt = long_exts->find(lpivot);
    assert( lxt != long_exts->end() ) ;
    
    Extender lnn = lxt->second;
    std::multimap<int, PathId> lposs = util::sortByValue<PathId,int>(lnn.positions);
    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "#Members:" << lnn.members.size() << "\n";
        for ( PathIdSet::iterator pt = lnn.members.begin(); 
              pt != lnn.members.end(); ++pt ) 
            std::cout << "l" << *pt << ":" << lnn.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << lnn.aligns.size() << "\n";
    }
    std::multimap<int, PathId>::iterator lit;
    for ( lit = lposs.begin(); lit != lposs.end(); ++lit ) {
        PathId lmember = lit->second;
        int    lpos = lit->first;
        assert( lpos >= 0 );
        AlignSummary lsum = findAlignSummary(lmember, lnn.aligns);
        if ( Param::verbose > 1 ) lsum.print(std::cout);
        if ( lsum.ilist.size() || lsum.dlist.size() ) {
            ct_gappy_extend++;
            gappy = true;
        }
        if ( Param::verbose > 1 && gappy ) std::cout << "gappy\n";


        assert( long_idmap->find(lit->second) != long_idmap->end() );

        PathId ppid = (*long_idmap)[lit->second];
        int    ppos = lit->first;
        if ( Param::verbose > 1 || Param::debug_flag ) {
            std::cout << "Long:" << lit->second << "\tPath:p" << ppid << "\tpos:" << ppos << "\n";
            std::cout << (*long_exts)[lit->second].sequence << "\n";
        }
        
        PathAligner aln;
        alignPath( aln, ppid );//, lsum );
        aln.adjust( ppos );
        laln.merge(aln);
    }
    if ( Param::verbose > 1 ) {
        std::cout << "Long extension placement\n";
        laln.printPlacement(std::cout);        
    }
}

void Placer::alignShortExtensionCluster( PathId ext_pid,
                                         bool &gappy,
                                         PathAligner &paln )
{
    ExtenderMap::iterator mt = tiny_exts->find(ext_pid);
    assert( mt != tiny_exts->end() );
    Extender ex = mt->second;

    if ( Param::pair_flag && ( Param::verbose > 1 || Param::debug_flag ) ) {
        std::cout << "\nShort overlapping path:l" << ext_pid << "\n";
        std::cout << ex.sequence << "\n";
        std::cout << "#Members:" << ex.members.size() << "\t";
        for ( PathIdSet::iterator pt = ex.members.begin(); 
              pt != ex.members.end(); ++pt ) 
            std::cout << "l" << *pt << ":" << ex.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << ex.aligns.size() << "\n";
    }
    
    std::multimap<int, PathId> posmap = util::sortByValue<PathId,int>(ex.positions);
    
    std::multimap<int, PathId>::iterator it;
    for ( it = posmap.begin();
          it != posmap.end(); ++it ) {
        PathId pid = it->second;
        int    pos = it->first;
        AlignSummary sum = findAlignSummary(pid, ex.aligns);
        if ( Param::verbose > 1 || Param::debug_flag ) std::cout << "Member:l" << pid << "\tpos:" << pos << "\n";

        assert( tiny_idmap->find(pid) != tiny_idmap->end() );
        PathId cid = (*tiny_idmap)[pid];
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Long-pid:c" << cid << "\n";
        
        PathAligner caln;
        alignLongExtensionCluster(cid, gappy, caln);
        if ( sum.dlist.size() || sum.ilist.size() ) {
            ct_gappy_extend++;
            gappy = true;
        }

        caln.adjust(pos);        
        paln.merge(caln);
    }
}

void Placer::alignLongExtensionCluster( PathId ext_pid,
                                     bool &gappy,
                                     PathAligner &paln )
{
    ExtenderMap::iterator mt = long_exts->find(ext_pid);
    assert( mt != long_exts->end() );
    Extender ex = mt->second;

    if ( Param::pair_flag && ( Param::verbose > 1 || Param::debug_flag ) ) {
        std::cout << "\nLong overlapping path:l" << ext_pid << "\n";
        std::cout << ex.sequence << "\n";
        std::cout << "#Members:" << ex.members.size() << "\t";
        for ( PathIdSet::iterator pt = ex.members.begin(); 
              pt != ex.members.end(); ++pt ) 
            std::cout << "l" << *pt << ":" << ex.positions[*pt] << " ";
        std::cout << "\n";
        std::cout << "#Aligns:" << ex.aligns.size() << "\n";
    }

    std::multimap<int, PathId> posmap = util::sortByValue<PathId,int>(ex.positions);
    
    std::multimap<int, PathId>::iterator it;
    for ( it = posmap.begin();
          it != posmap.end(); ++it ) {
        PathId pid = it->second;
        int    pos = it->first;
        AlignSummary sum = findAlignSummary(pid, ex.aligns);
        if ( Param::verbose > 1 || Param::debug_flag ) std::cout << "Member:l" << pid << "\tpos:" << pos << "\n";

        assert( long_idmap->find(pid) != long_idmap->end() );
        PathId cid = (*long_idmap)[pid];
        if ( Param::verbose > 1 || Param::debug_flag ) 
            std::cout << "Cluster-pid:c" << cid << "\n";
        
        PathAligner caln;
        alignCluster(caln, cid, pos, sum, gappy);

        if ( sum.dlist.size() || sum.ilist.size() ) {
            ct_gappy_extend++;
            gappy = true;
        }

        caln.adjust(pos);        
        paln.merge(caln);
    }
}

void Placer::alignCluster( PathAligner &caln, 
                           PathId pid, 
                           int cpos,
                           AlignSummary &sum, 
                           bool &gappy )
{
    AlignPosList ins = sum.ilist;
    AlignPosList del = sum.dlist;

    Clusters::iterator it = clusters->find(pid);
    assert( it != clusters->end());
    
    assert( cluster_idmap->find(pid) != cluster_idmap->end() );
    PathId opid = (*cluster_idmap)[pid];
    
    if ( Param::verbose > 1 || Param::debug_flag ) 
        std::cout << "Center:c" << pid << "\tPath:" << opid << "\n";

    
    //alignPath( caln, pid, sum );

    IntPair stops = (*read_aligns)[opid].getStops();
    IntPair trims = (*read_aligns)[opid].getTrims();

    if ( Param::verbose > 1 || Param::debug_flag ) 
        std::cout << "pid:c" << pid << " stops:" << stops.first << stops.second
                  << " trims:" << trims.first << trims.second << "\n";

    PathIdList::iterator pt;
    AlignSummaryList::iterator at;
    
    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "#Members:" << it->second.members.size() << "\n";
        PathIdList::iterator pt;
        AlignSummaryList::iterator at;
        for ( pt  = it->second.members.begin(), at =  it->second.summarys.begin();
              pt != it->second.members.end(),   at != it->second.summarys.end();
              ++pt, ++at ) {
            std::cout << "c" << *pt << ":" << at->lgap.second << " ";
        }
        std::cout << "\n";
    }
    
    for ( pt  = it->second.members.begin(), at =  it->second.summarys.begin();
          pt != it->second.members.end(),   at != it->second.summarys.end();
          ++pt, ++at ) {
        PathId id = *pt;
        AlignSummary as = *at;

        
        assert( cluster_idmap->find(id) != cluster_idmap->end() );
        PathId oid = (*cluster_idmap)[id];

        if ( Param::verbose > 1 || Param::debug_flag ) {
            std::cout << "Path:p" << id << "\torig:" << oid << "\t";
            std::cout << (*read_aligns)[oid].getReference() << "\n";
        }
        
        PathAligner aln;
        alignPath( aln, oid );

        int ppos = as.lgap.second;
        if ( Param::verbose > 1 || Param::debug_flag ) 
            printf("shift-path:%d\n", ppos);
        aln.adjust( ppos );
        
        caln.merge(aln);

        if ( as.dlist.size() || as.ilist.size() ) {
            ct_gappy_merge++;
            gappy = true;
        }
    }
}

void Placer::alignPath( PathAligner& paln,
                        PathId pid )
{
    // if ( Param::verbose || Param::debug_flag ) {
    //     std::cout << "range:" << sum.range.first << " " << sum.range.second << "\n";
    //     std::cout << "lgap:"  << sum.lgap.first << " " << sum.lgap.second << "\n";
    //     std::cout << "egap:"  << sum.egap.first << " " << sum.egap.second << "\n";
    //     std::cout << "dlist:" << sum.dlist.size() << "\n";
    //     std::cout << "ilist:" << sum.ilist.size() << "\n";
    // }

    PathType     path;
	KmerArray    kmers;
	SuffixBounds bounds;

    assert( (*read_aligns).find(pid) != (*read_aligns).end() );
    ReadAligner raln = (*read_aligns)[pid];
    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "PathId:p" << pid << "\n";
        std::cout << raln.getReference() << "\n";
        std::cout << "stops:" << raln.getLStop() << " " << raln.getRStop() << "\n";
        std::cout << "trims:" << raln.getLTrim() << " " << raln.getRTrim() << "\n";
    }
    paln = PathAligner(raln);
}

bool Placer::findPairRead( ReadId rid, PathId pid )
{
    PathToNumMap::iterator pt = Path2Ids.find(pid);
    assert( pt != Path2Ids.end() );
 
   ReadPlacementList *places = Alns[pt->second].getPlacements();
   ReadPlacementList::iterator it;
   for ( it = places->begin(); it != places->end(); ++it )
       if ( it->rid == rid ) return true;

   return false;
}

void Placer::check()
{
    std::cout << "\nPath sanity check ...\n";
    int bad = findBadPaths();
    std::cout << "#invalid paths:" << bad << "\n";    
    
    //countInvalidPaths();
}

int Placer::findBadPaths()
{
    progress = Progress( 1.0, 0, npaths, mytime() );

    std::vector<size_t> bad_counts(4,0);

    size_t i;
    int why;
#pragma omp parallel for schedule(dynamic, 1) if(Param::ncpus>1) private(i) num_threads(Param::ncpus)    
    for ( i = 0; i < npaths; i++ ) {
        if ( Param::debug_id != -1 && Param::debug_id != Id2Paths[i] ) continue;
        if ( !validPath(i, why) ) {
            bad_paths[i] = true;
            if ( Param::verbose > 1 ) {
                std::cout << "Invalid path:\tId:" << i 
                          << "\tpid:" << Id2Paths[i] << "\t"
                          << why << "\n";
                
#pragma omp critical
                {
                    bad_counts[why]++;
                }
            }
        }
#pragma omp critical
        {
            double prev = progress.ratio;
            progress.count++;
            progress.showProgress();

            if ( ( Param::verbose >= 1 ) && progress.ratio > prev )
                countInvalidPaths();
        }
        
    }
    
    int num = 0;
    for ( size_t i = 0; i < npaths; i++ ) 
        if ( bad_paths[i] ) num++;

    return num;
}

void Placer::countInvalidPaths()
{

    size_t ct_msa_error = ct_msa_empty + ct_msa_lerr + ct_msa_rerr + ct_msa_hole;
    size_t ct_con_error = ct_con_lerr + ct_con_rerr + ct_con_hole;
    size_t ct_bad = ct_empty_ref + ct_msa_error + ct_con_error;

    std::cout << "Reference sequence trimmed:" << "(left:" << ct_msa_ltrim << ")\n";
    std::cout << "Consenus sequence trimmed:" << "(left:" << ct_con_ltrim << "," << "right:" << ct_con_rtrim << ")\n";

    if ( ct_bad == 0 ) return;
    std::cout << "Invalid paths:" << ct_bad << "\n";
    if ( ct_empty_ref > 0 ) 
        std::cout << "- Empty reference:" << ct_empty_ref << "\n";
    if ( ct_msa_error > 0 ) {
        std::cout << "- MSA errors:" << ct_msa_error << "\t";
        std::cout << "(empty sequence:" << ct_msa_empty << ", ";
        std::cout << "leading sequence error:" << ct_msa_lerr << ", ";
        std::cout << "trailing sequence error:" << ct_msa_rerr << ", ";
        std::cout << "coverage hole:" << ct_msa_hole << ")\n";
    } 
    if ( ct_con_error > 0 ) {
        std::cout << "- Consensus errors:" << ct_con_error << "\t";
        std::cout << "(leading sequence error:" << ct_con_lerr << ", ";
        std::cout << "trailing sequence error:" << ct_con_rerr << ", ";
        std::cout << "coverage hole:" << ct_con_hole << ")\n";
    }
}

bool Placer::validPath(int i, int &why)
{
    PathId pid = Id2Paths[i];
    std::string seq = Alns[i].getSequence();
    ReadPlacementList *places = Alns[i].getPlacements();

    if ( Param::verbose > 1 || Param::debug_flag ) {
        std::cout << "\nID:"<< i << "\n";
        std::cout << "pid:" << pid << "\n";
        std::cout << "seq:" << seq << "\n";
        std::cout << "len:" << seq.size() << "\n";
        std::cout << "# reads:" << places->size() << "\n";
    }
    if ( seq.size() == 0 ) {
        if ( Param::verbose > 1 || Param::debug_flag ) std::cout << "Empty sequence\n";
        why = EMPTY_SEQUENCE;
        ct_empty_ref++;
        return false;
    }

    if ( places->size() == 0 ) {
        if ( Param::verbose >= 1 ) {
            std::cout << "\nID:"<< i << "\n";
            std::cout << "pid:" << pid << "\n";
            std::cout << "seq:" << seq << "\n";
            std::cout << "len:" << seq.size() << "\n";
            std::cout << "# reads:" << places->size() << "\n";
            std::cout << "Zero reads\n";
        }
        why = ZERO_READS;
        ct_zero_reads++;
        return false;
    }

    // bool valid = true;   
    // for ( ReadPlacementList::iterator it = places->begin(); it != places->end(); ++it ) {
    //     if ( it->read_pos < 0 || it->read_pos >= strlen(seqs[it->rid]) ) {
    //         std::cout << "Invalid rpos\tread:" << it->rid << "\t" << seqs[it->rid] << "\trpos:" << it->read_pos << "\tppos:" << it->path_pos << "\n";
    //         valid = false;
    //     }
    //     if ( it->path_pos < 0 || it->path_pos >= seq.size() ) {
    //         std::cout << "Invalid ppos\tread:" << it->rid << "\t" << seqs[it->rid] << "\trpos:" << it->read_pos << "\tppos:" << it->path_pos << "\n";
    //         valid = false; 
    //     }
    //     std::list<Mismatch>::iterator mt;
    //     for ( mt = it->ilist.begin(); mt != it->ilist.end(); ++mt ) {
    //         if ( mt->qry_pos < 0 || mt->qry_pos >= strlen(seqs[it->rid]) ) {
    //             std::cout << "Invalid INS\tread:" << it->rid << "\t" << seqs[it->rid] << "\trpos:" << mt->qry_pos << "\tppos:" << mt->ref_pos << "\n";
    //             valid = false;
    //         }
    //         if ( mt->ref_pos < 0 || mt->ref_pos >= seq.size() ) {
    //             std::cout << "Invalid INS\tread:" << it->rid << "\t" << seqs[it->rid] << "\trpos:" << mt->qry_pos << "\tppos:" << mt->ref_pos << "\n";
    //             valid = false;
    //         }
    //     }
    //     for ( mt = it->dlist.begin(); mt != it->dlist.end(); ++mt ) {
    //         if ( mt->qry_pos < 0 || mt->qry_pos >= strlen(seqs[it->rid]) ) {
    //             std::cout << "Invalid DEL\tread:" << it->rid << "\t" << seqs[it->rid] << "\trpos:" << mt->qry_pos << "\tppos:" << mt->ref_pos << "\n";
    //             valid = false;
    //         }
    //         if ( mt->ref_pos < 0 || mt->ref_pos >= seq.size() ) {
    //             std::cout << "Invalid DEL\tread:" << it->rid << "\t" << seqs[it->rid] << "\trpos:" << mt->qry_pos << "\tppos:" << mt->ref_pos << "\n";
    //             valid = false;
    //         }
    //     }
    // }

    // if ( valid == false ) {
    //     std::cout << "\nID:"<< i << "\n";
    //     std::cout << "pid:" << pid << "\n";
    //     std::cout << seq << "\n";
    //     std::cout << "len:" << seq.size() << "\n";
    //     std::cout << "# reads:" << places->size() << "\n";
    //     std::cout << "Invalid position\n";
    // }


    MSA msa( &Alns[i], seqs) ;

    // Do not generate consensus here.
    //Profile pro(&msa, AA_TYPE);
    std::string ref = Alns[i].getSequence();
    //std::string con = pro.makeConsensus();
    std::string pad = msa.getPivot();
    if ( pad.size() == 0 ) {
        why = EMPTY_SEQUENCE;
        ct_msa_empty++;
        return false;
    }

    // trim
    int s = 0, e = (int)pad.size()-1;
    while ( pad[s] == '.' ) s++;
    if ( s >= e ) {
        why = TRIM_ERROR;
        ct_msa_lerr++;
        return false;
    }
    while ( pad[e] == '.' ) e--;
    if ( e <= s ) {
        why = TRIM_ERROR;
        ct_msa_rerr++;
        return false;
    }

    if ( Param::verbose >= 1 || Param::debug_flag ) {
        if ( s > 0 ) 
            std::cout << "Seqence trimmed\tId:" << i 
                      << "\tpid:" << pid 
                      << "\tltrim:" << s << "\n";

        // if ( e < (int)pad.size()-1 ) 
        //     std::cout << "Seqence trimmed\tId:" << i 
        //               << "\tpid:" << pid 
        //               << "\trtrim:" << (int)pad.size()-1-e << "\n";
    }
    
    
    // find hole
    for ( int p = s; p<=e; p++ ) {
        if ( pad[p] == '.' ) {
            why = COVERAGE_HOLE;
            ct_msa_hole++;
            return false;
        }
    }

    // NO adjust in right side
    /*
 ID:751
        Path:117788
        0          .    :    .    :    .    : 
 reference      VSTKPQPAAATEDFLTDCTVHFPRIPF....
   7816746      VSTKPQPAAATEDFLTDCTVHFPRIPFGMLF
   5070296      VSTKPQPAAATEDFLTDCTVHFPRIPFG...
   3246595      VSTKPQPAAATEDFLTDCTVHFPRIPFG...
   2369212      VSTKPQPAAATEDFLTDCTVHFPR.......
  16199822      VSTKPQPAAATEDFLTDCTVHFPR.......
    */

    // // adjust sequence
    // if ( e < (int)pad.size()-1 ) {
    //     pad = pad.substr(0, e);
    //     if ( e < (int)ref.size() ) {
    //         std::string old = ref;
    //         ref = ref.substr(0, e);
    //         if ( Param::verbose ) {
    //             std::cout << "old ref:" << old << "\n";
    //             std::cout << "new ref:" << ref << "\n";
    //         }
    //     }

    //     int d = (int)pad.size()-1 - e;
    //     Alns[i].setRTrim( d+Alns[i].getRTrim() );
    // }

    // adjust alignment
    if ( s > 0 ) {
        ct_msa_ltrim++;
        std::string old = ref;
        Alns[i].adjust(-1*s);
        pad = pad.substr(s);
        ref = ref.substr(s);
        if ( Param::verbose > 1 || Param::debug_flag ) {
            std::cout << "old ref:" << old << "\n";
            std::cout << "new ref:" << ref << "\n";
        }
        
        Alns[i].setLTrim( s+Alns[i].getLTrim() );
    }

    // CAUTION: No not remove gap.
    // It makes alignment awkward if realign
    // for ( int p = (int)pad.size()-1; p >= 0; p-- ) 
    //     if ( pad[p] == '-' ) con.erase(p,1);

    //if ( Param::verbose || Param::debug_flag ) std::cout << "Consensus:" << pad << "\n";

    if ( ref != Alns[i].getSequence() ) 
        Alns[i].setSequence(ref); 
    //Alns[i].setConsensus(pad);

    Profile pro(&msa, AA_TYPE);
    std::string con = pro.makeConsensus();

    /*
    Mon Feb  2 14:14:08 EST 2015
    Trim off leading and trailing dots
    This happens when an augmented path is re-aligned.
    Some of reads may be dropped because of low sequence simlarity
    against new path.
    */
    s = 0, e = (int)con.size()-1;
    while ( con[s] == '.' ) s++;
    if ( s >= e ) {
        why = TRIM_ERROR;
        ct_con_lerr++;
        return false;
    }
    while ( con[e] == '.' ) e--;
    if ( e <= s ) {
        why = TRIM_ERROR;
        ct_con_rerr++;
        return false;
    }

    // find hole
    for ( int p = s; p<=e; p++ ) {
        if ( con[p] == '.' ) {
            why = COVERAGE_HOLE;
            ct_con_hole++;
            return false;
        }
    }

    if ( e < (int)con.size()-1 ) {
        ct_con_rtrim++;
        int m = (int)con.size()-1 - e;
        con = con.substr(0, e+1);
        assert(m < (int)ref.size()-1);
        int r = (int)ref.size()-1 - m;
        ref = ref.substr(0,r);

        Alns[i].setRTrim( m+Alns[i].getRTrim() );    
        //Alns[i].dropGaps(e);
        
    }
    if ( s > 0 ) {
        ct_con_ltrim++;
        assert(s<(int)con.size());
        con = con.substr(s);
        assert(s<(int)ref.size());
        ref = ref.substr(s);
        Alns[i].setLTrim( s+Alns[i].getLTrim() );    
        Alns[i].adjust(-1*s);
    }

    // Update reference sequence as well if consensus is changed
    if ( ref != Alns[i].getSequence() ) 
        Alns[i].setSequence(ref);

    size_t csize = con.size(), rsize = ref.size();
    assert( csize > 0 && rsize > 0 );
    if ( csize != rsize && Param::verbose > 1 ) {
        std::cout << "Path:" << i << "\tLength diff\t";
        std::cout << "Consensus:" << csize << "\tReference:" << rsize << "\n";
        std::cout << "Consensus:" << con << "\n";
        std::cout << "Reference:" << ref << "\n";
        std::cout << "Original :" << seq << "\n";
    }
    size_t end = std::min( csize, rsize ) - 1;
    bool found = Alns[i].dropGaps(end); // gaps will be trimmed
    if ( found && Param::verbose > 1 ) {
        std::cout << "Path:" << i << "\tFound out of bound gaps\n";
        std::cout << "Consensus:" << con << "\n";
        std::cout << "Reference:" << ref << "\n";
        std::cout << "Original :" << seq << "\n";
    }
    Alns[i].setConsensus(con);
    Alns[i].setDirty(false);
    if ( pro.replacedStopCodon() ) {
        if ( Param::verbose > 1 ) 
            std::cout << "Id:" << i 
                      << "\tpid:" << Id2Paths[i] 
                      << "\tstop codon replaced\n";
    }
    
    return true;
}
