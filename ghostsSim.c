#include <stdio.h>
#include <stdlib.h>
#include <time.h>



//simulation type switches
#define simulation_type_default 0 //default: constant birth rate


//rate parameters
#define b_default 1 //birth
#define d_default 0.5 // death 
#define m_default 0.01 //out-migration .01 -> .1
#define q_default 0.5 //frequency of immigrant altruists (0,1) by .1

//simulation constants
#define n0_default 500 //initial population size 
#define a0_default 250 //initial number of altruists size
#define K_default 1050 //carring capacity
#define M_default 1000 //size at which populations fission  500 -> 4k
#define T_default 100000 //total number of turns for which the simulation should run
#define N_default 1000 // Total number of populations we need to track


//initialise variables to defaults

//switches
int simulation_type = simulation_type_default;


//rates
double b = b_default; //birth rate
double d = d_default; //death rate
double m = m_default; //migration rate
double q = q_default; //frequency of ghostly altruists
double phi = 0; //in-migration rate; we set it dynamically later.

//ints
int N = N_default;
int M = M_default;
int T = T_default;
long K = K_default;
	//int G = T_default;


//initial sizes
int a0 = a0_default;
int n0 = n0_default;


//temporary variables
int Q;  // integer that must be returned by rand() for U(0,1) > q
long tempA; //temporary counter for number of altruists while fissioning
long tempN; //temporary counter for number of non-altruists while fissioning
long parent; // a randomly selected individual
double r; //random number



//simulation tracking variables 
long *ns; //vector of population sizes
long *as; //vector of altruists in each population




//simulation functions
void simulate();
void event(long *n, long *a);
void reproduce(long *n, long *a);
void expire(long *n, long *a);
void migrateIN(long *n, long *a);
void fission(long *n, long *a);
void printResults();
void usage();


//main loop
int main(int argc, char *argv[]){ 
	srand(time(NULL)); //crank up the random number generator
	
	int i;
	
	for (i = 1; i < argc; i++) { //parse command line arguments    
		if (argv[i][0] == '-') {
				 if (argv[i][1] == 'h') usage(argv[0]);
			else if (argv[i][1] == 'q') q = atof(argv[i+1]);
			else if (argv[i][1] == 'b') b = atof(argv[i+1]);
			else if (argv[i][1] == 'd') d = atof(argv[i+1]);
			else if (argv[i][1] == 'm') m = atof(argv[i+1]);
			else if (argv[i][1] == 'i') phi = atoi(argv[i+1]);
			else if (argv[i][1] == 'n') n0 = atoi(argv[i+1]);
			else if (argv[i][1] == 'a') a0 = atoi(argv[i+1]);
			else if (argv[i][1] == 'M') M = atoi(argv[i+1]);
			else if (argv[i][1] == 'T') T = atoi(argv[i+1]);
			else if (argv[i][1] == 'N') N = atoi(argv[i+1]);
			else if (argv[i][1] == 'K') K = atoi(argv[i+1]);
			else if (argv[i][1] == 'S') simulation_type = atoi(argv[i+1]);
		}
	}
	
	if(phi <= 0){ //in-migration wasn't set, we set it dynamically here
		phi = m * M * 0.75;
	}
	
	ns = malloc(N * sizeof(long));
	as = malloc(N * sizeof(long));
	
	Q = (int)(q*RAND_MAX); //integer that must be returned by rand() for real > q
	
	//initialize N populations
	for (i = 0; i < N; i++){ 
		ns[i] = n0;	
		as[i] = a0;
	}
		
	simulate(); //run simulation	
	printResults(); //output results to stdout
	exit(0);
}


void usage(char *name){
	fprintf (stderr, "\
	Usage: %s [OPTIONS]\n \
	Simulates a Polya's Urn model of population growth, with population fissioning (each individual leaves with pr = 0.5), and the following parameters: \n\
	\n\
	Transition Rates \n\
	-b   Birth rate (b * n; default: 1) \n\
	-d   Death rate (d * n; default: 0.5) \n\
	-m   out-migration rate (m * n; default: 0.01) \n\
	-i   Absolute in-migration rate (phi) (default: m*M*.75) \n\
	-q   In-migrating altruist proporiton (default: 0.5) \n\
	\n\
	Constants \n\
	-M   Size at which population fissions (default: 1000) \n\
	-K   Carrying capacity, for dynamic birth or death rates (default: 1050) \n\
	-T   Number of events (turns) to simulate (default: 100000) \n\
	-N   Number of populations to simulate (default: 1000) \n\
	\n\
	Simulation variants \n\
	-S   Sim. type:\n\
	         0: static birth and death rates (default)\n\
	         1: frequency dependent death rate: b*n/k\n\
	         2: frequency dependent birth rate: 1-((1-d)*n/k)\n\
	\n\
	Initialisation \n\
	-n   Initial size of each population (default: 500) \n\
	-a   Initial number of altruists in each population (default: 250) \n\
	\n\
	Example: %s -t 1 -d 1 -m 0.05 -q .2 -N 100 -T 1000000\
	", name, name);
	exit(0);
}


void migrateIN(long *n, long *a){ //someone migrates in, it's an altruist with probability q
	if(rand() <= Q){ (*a)++; }
	(*n)++;
}

void reproduce(long *n, long *a){ //determines who is born and calls birth(int altruist)	
	parent = 1+(rand() % *n);
	if (parent <= *a) (*a)++; //an alrtuist breeds			
	(*n)++; //someone breeds
}

void expire(long *n, long *a){ //kills one random individual in the population	
	parent = 1+(rand() % *n);
	if (parent <= *a) (*a)--; //an altruist dies
	(*n)--;
}

void fission(long *n, long *a){ //kills one random individual in the population	
	if(*n > M){
		tempA = tempN = 0;
		int i; //variable for iterating
		for (i = 0; i < *a; i++){
			tempA += (rand() %2);
		}	
		for (i = *a; i < *n; i++){
			tempN += (rand() %2);
		} 
		*a = tempA;
		*n = tempA+tempN;
	}
}


void event(long *n, long *a){ //one turn of the simulation		
	r = ((double) rand() / RAND_MAX); //random number between 0 and 1	
	r *=  ((*n) * (b+d+m) + phi);  //random number between 0 and n(b + d + m) + phi
	if(r<=phi) migrateIN(n, a); //an in-migrantion
	else{
		if ( r <= ( b * (*n) + phi  ) ) reproduce(n, a); //a birth
		else expire(n, a); //a death or out-migration
	}
	//printf("--- %ld, %ld, %lf, %lf ---", *n,k, d, bdm);
	//if( t % 100000 == 0 ) fprintf(stderr, "%d\n", t); //uncomment to print progress 
}

void simulate(){//We select a loop one based on the simulation type parameters, to avoid having to recheck them on each event.
	int i; //variable for iterating
	int t; //the current turn of the simulation
	if (simulation_type == 2) { //Frequency dependent births, constant deaths; implies constant d
		for (t = 0; t < T; t++){
			for (i = 0; i < N; i++){
				fission(&ns[i], &as[i]);
				b = 1 - (ns[i] /  (double)K) * (1-d);
				event(&ns[i], &as[i]);
			}
		}
	} else if (simulation_type == 1) { //Frequency dependent deaths, constant births	
		for (t = 0; t < T; t++){
			for (i = 0; i < N; i++){
				fission(&ns[i], &as[i]);
				d = b * (ns[i] /  (double)K);
				event(&ns[i], &as[i]);
			}
		}
	} else { //constant death and birth, default
		for (t = 0; t < T; t++){
				for (i = 0; i < N; i++){
					fission(&ns[i], &as[i]);
					event(&ns[i], &as[i]);				
			}
		}
	}
}


void printResults(){ //prints simulation varaibles to stdout
	char *format = "%d, %d, %f\n";
	printf("%s, %s, %s\n", "n","a","p");
	int i;
	for (i = 0; i < N; i++) printf(format, ns[i],as[i], ((double)as[i]/ns[i]) );
}
