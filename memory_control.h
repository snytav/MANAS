/*
 * memory_control.h
 *
 *  Created on: Nov 16, 2016
 *      Author: snytav
 */

#ifndef MEMORY_CONTROL_H_
#define MEMORY_CONTROL_H_


int get_mem_used_sizeMb();
int get_mem_used_various(int nt,char *where);
int memoryAllocationLog(size_t size,char *fname,int line_number);


#endif /* MEMORY_CONTROL_H_ */
