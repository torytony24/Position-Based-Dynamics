#pragma once
#include <string>
#include <windows.h>
#include <winbase.h>

class Performance {
public:
	LARGE_INTEGER freq, start, end;

	std::string name;
	Performance() {};
	Performance(std::string _name) {
		name = _name;
	}

	inline void Start() {
		QueryPerformanceCounter(&start);
		QueryPerformanceFrequency(&freq);
	}

	inline double End() {
		QueryPerformanceCounter(&end);
		double time = ((double)(long double)(end.QuadPart - start.QuadPart) / (long double)freq.QuadPart);

		double msTime = 1000 * time;
		//cout << name << "  : " << msTime << "ms" << endl;
		printf("%-30s : %10.3f ms\n", name.c_str(), msTime);
		return msTime;
	}
};
