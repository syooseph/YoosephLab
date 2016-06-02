#include "sectioner.h"

Sectioner::Sectioner()
{
    section.count = 0;
}

Sectioner::Sectioner( std::string *q, std::string *s )
{
    init( q, s );
}

void Sectioner::init( std::string *q,
                      std::string *s )
{
    query = q;
    sbjct = s;

    section = Section(0,0,0,0,0);
}

