#include "stdafx.h"
#include "MemorySurveillance.h"

using namespace cwg;
using namespace std;
//------------------------------------------------------------------------------
bool MemorySurveillance::AdjustPurview() const
{
	TOKEN_PRIVILEGES TokenPrivileges;
	BOOL bRet;
	HANDLE hToken;

	LookupPrivilegeValue(NULL, SE_DEBUG_NAME, &TokenPrivileges.Privileges[0].Luid);   
	OpenProcessToken(GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES, &hToken);

	TokenPrivileges.PrivilegeCount = 1;   
	TokenPrivileges.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;

	bRet = AdjustTokenPrivileges(hToken, FALSE, &TokenPrivileges, 0, NULL, NULL);

	CloseHandle(hToken);
	return bRet == TRUE;
}
//------------------------------------------------------------------------------
PROCESS_MEMORY_COUNTERS MemorySurveillance::GetMemoryUsage() const
{
	HANDLE handle=GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(handle,&pmc,sizeof(pmc));
	return pmc;
}
//------------------------------------------------------------------------------
MemorySurveillance& MemorySurveillance::GetInstance()
{
	static MemorySurveillance kInstance;
	return kInstance;
}
//------------------------------------------------------------------------------
void cwg::ShowPagefileMemoryUsage(std::ostream& out)
{

	double dMemoryUsage = MemorySurveillanceInstance.GetMemoryUsage().PagefileUsage/1024.0/1024.0;
	out << "Pagefile Usage: " << dMemoryUsage << " MB" << endl;
}
//------------------------------------------------------------------------------