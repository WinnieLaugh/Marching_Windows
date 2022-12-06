#ifndef _MEMORYSURVEILLANCE_H_
#define _MEMORYSURVEILLANCE_H_

#include <stdio.h>    
#include <windows.h> 
#include <Psapi.h>
#include <iostream>
#include <string>
#include <ostream>

#pragma comment(lib, "psapi.lib")

namespace cwg
{
	class MemorySurveillance
	{
	public:
		bool							AdjustPurview() const;
		PROCESS_MEMORY_COUNTERS			GetMemoryUsage() const;
		static MemorySurveillance&		GetInstance();
	
	private:
		MemorySurveillance()			{ AdjustPurview(); };
		MemorySurveillance( const MemorySurveillance& ){};
		MemorySurveillance& operator = ( const MemorySurveillance& ){};
	};

#define MemorySurveillanceInstance MemorySurveillance::GetInstance()

	void								ShowPagefileMemoryUsage(std::ostream& out = std::cout);

}

#endif