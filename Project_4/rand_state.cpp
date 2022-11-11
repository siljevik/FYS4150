// member-variable and -function defined here 
#include "MCMC_spin.hpp" 

// In this file we generate a random integer (state) -1 or +1

// hentet fra eksemplel
// https://github.com/anderkve/FYS3150/blob/master/code_examples/random_number_generation/src/omp_rng.cpp

// Initialize the omp_rng_container according to the number of threads
  void initialize_omp_rng_container(unsigned int base_seed)
  {

    // Let thread 0 add one generator to omp_rng_container for each thread
    #pragma omp parallel
    {
      // Which thread am I?
      int thread_id = omp_get_thread_num();    

      // Let thread 0 do the admin work
      if (thread_id==0)
      {
        // How many threads available in total?
        int n_threads = omp_get_num_threads();

        // If the base_seed isn't set, use the system clock
        if (base_seed == -1)
        {
          base_seed = chrono::system_clock::now().time_since_epoch().count();
        }

        // Add one generator per thread. Seed each generator with base_seed + thread number
        for (int i = 0; i < n_threads; i++)
        {
          mt19937 my_generator(base_seed + i);
          omp_rng_container.push_back(my_generator);
        }

      }

      // Don't let any thread run past this point and try using 
      // the random number generators until thread 0 is done creating them
      #pragma omp barrier

    } // End parallel region

  }
  

int get_random_int_m1_p1(){
  static uniform_int_distribution<int>distribution(-1,1);

  // sample a number using the generator for the current thread
  int thread_id = omp_get_thread_num();
  int r = distribution(omp_rng_container.at(thread_id));

  return r;
};