/*  Link in this file for random number generation with rand()
     and seeding from the clock  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

	double
ran1()
{
	int rand();
	return( rand()/(RAND_MAX+1.0)  );
}


	void seedit( const char *flag)
{
	// FILE *fopen(), *pfseed;
	unsigned int seed2 ;

    if( flag[0] == 's' ) {
        seed2 = time(NULL);
        srand(seed2);
        printf("\n%d\n", seed2 );    
    }

}

	int
commandlineseed( char **seeds)
{
	unsigned int seed2 ;
    // void srand(unsigned int seed);

	seed2 = atoi( seeds[0] );

	// printf("\n%d\n", seed2 );    

	srand(seed2); 
	return(1);
}

