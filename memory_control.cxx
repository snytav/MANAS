/*
 * memory_control.cxx
 *
 *  Created on: Nov 16, 2016
 *      Author: snytav
 */


#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi_shortcut.h"
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sys/types.h>

#include <string>
#include <sstream>
#include <vector>

using namespace std;

unsigned long int total_memory_volume = 0;




int GetDeviceMemory(size_t *m_free,size_t *m_total);

int get_mem_used_sizeMb()
{
  int i = 0;
  struct rusage r_usage;
  getrusage(RUSAGE_SELF,&r_usage);

  return (r_usage.ru_maxrss/1024/1024);
}

typedef struct {
    unsigned long size,resident,share,text,lib,data,dt;
} statm_t;

void read_off_memory_status(statm_t& result)
{
  unsigned long dummy;
  const char* statm_path = "/proc/self/statm";

  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
    &result.size,&result.resident,&result.share,&result.text,&result.lib,&result.data,&result.dt))
  {
    perror(statm_path);
    abort();
  }
  fclose(f);
}

void read_free_total_memory(unsigned long &m_total, unsigned long &m_free)
{
  unsigned long dummy;
  const char* statm_path = "/proc/meminfo";
  char str[100];

  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  if(fgets(str,100,f)  != NULL)
  {
      m_total = atoi(str + 10);
  }
  if(fgets(str,100,f)  != NULL)
  {
      m_free = atoi(str + 10);
  }
  fclose(f);
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    int n = 0;
    while (std::getline(ss, item, delim)) 
    {
        if(item != "")
	{
           elems.push_back(item);
	   n++;
	}
    }
}

float get_fraction_memory_used()
{
        char s[100],fname[100],str[100];
	FILE *f;
	string str_copy,frac;
	vector<string> vs;
	float frac_used = 0.0;
	
        // get fraction of the memory used by the program through ps -p <pid> v
	pid_t pid = getpid();
	sprintf(s,"ps -p %d v 2>&1 > mem_out_%d",pid,pid);
	system(s);
	
	sprintf(fname,"mem_out_%d",pid);
	
        if((f = fopen(fname,"rt")) == NULL) return -1.0;
	
	if(fgets(str,100,f)  != NULL)
	{
	    if(fgets(str,100,f)  != NULL)
	    {
	       str_copy =str;
	       split(str_copy,' ',vs);
	       
	    }
	}
	
	fclose(f);
	
	if(vs.size() > 8)
	{
	   frac = vs[8];
	   frac_used = atof(frac.c_str());
	}
	else
	{
	   frac_used = 0.0;
	}
	
	
	
	return frac_used;
}

int get_mem_used_various(int nt,char *where)
{
	statm_t result;
	FILE *f;
	char s[100];
	size_t m_free,m_total,m_free1,m_used,m_total1;
        struct sysinfo info;
	
	// get free memory through cat /proc/meminfo
	read_free_total_memory(m_total1,m_free1);
    
	
	
	float frac_used = get_fraction_memory_used();
	
	m_used = frac_used*m_total1*0.01;
	
    

	sprintf(s,"memory_usage_rank%05d.dat",getRank());

	if((f = fopen(s,"at")) == NULL) return 1;

	 sysinfo(&info);
	 GetDeviceMemory(&m_free,&m_total);

	read_off_memory_status(result);

	fprintf(f,"nt %8d %20s used %5.2f percent of %10ld Kb total = %10ld Kb free %10ld Kb (meminfo) RAM free %10lu USED: size %5ld resident %5ld share %5ld text %5ld lib %5ld data %5ld dt %5ld GPU total %10lu free %10lu (in Kb, from /proc/self/statm )\n",
			nt,
			where,
	                frac_used,m_total1,m_used,m_free1,
			info.freeram/1024,
			result.size/1024,
			result.resident/1024,
			result.share/1024,
			result.text/1024,
			result.lib/1024,
			result.data/1024,
			result.dt/1024,
			m_total/1024,
			m_free/1024
			);
	fclose(f);
	return 0;
}


int memoryAllocationLog(size_t size,string fname,int line_number)
{
    FILE *f;
    char s[100];
    int rank;
    size_t m_total1,m_free1,m_used;
    
    total_memory_volume += (unsigned long int)size;
    
//     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    sprintf(s,"memory_usage_rank%05d.dat",getRank());
    
    // get free memory through cat /proc/meminfo
	read_free_total_memory(m_total1,m_free1);
    
	
	
	float frac_used = get_fraction_memory_used();
	
	m_used = frac_used*m_total1*0.01;

    if((f = fopen(s,"at")) == NULL) return 1;
    
    fprintf(f,"                         %10lu at %20s line %5d total(Kb) %10ld from system %10ld \n",size,fname.c_str(),line_number,total_memory_volume/1024,m_used);
	
    fclose(f);
    
    return 0;
}


