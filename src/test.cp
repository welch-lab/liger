#include <math.h>
#include <iostream> 
#include <string>

char int2char(int digit)
{
	return digit+48;
}
void num2string_reverse(int num,std::string& str)
{
	str = "";
	while(num > 0)
	{
		str += (char)((num%10) + 48);
		num/=10;
	}
}

int intCharLength(int num)
{
	return floor(log10(num)+1);
}


int main(int argc,char* argv[])
{
	int num = 9348;
	std::cout << intCharLength(num) << std::endl;
	std::cout << num/10
}