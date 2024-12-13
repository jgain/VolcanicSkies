/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;*/
#include <iostream>
#include "sun.h"
#include "sunEngine.h"

typedef double valType;
using namespace PWM::Engine;
void testConstructors(){
	auto a = sunEngine<valType>(60.f, true);
	if (a.getSunList().size() == 0)
		std::cout << "\033[1;32mConstructor works as expected.\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! Constructor does not perform as expected!\033[0m" << std::endl;

}

void testMiscellany(){
	auto a = sunEngine<valType>(21600.f, true);
	std::shared_ptr<PWM::Model::sun<valType>> s = std::make_shared<PWM::Model::sun<valType>>(PWM::Model::sun<valType>("../resources/SunSol.json"));
	a.addSun(s);
	a.step();
	if (s->getApparentRightAscension() > 89. && s->getApparentRightAscension() < 91.)
		std::cout << "\033[1;32mstepApparentRightAscension() works as expected with a value of " << s->getApparentRightAscension() << ".\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! stepApparentRightAscension() does not perform as expected with a value of " << s->getApparentRightAscension() << "!\033[0m" << std::endl;
	
	for (int i = 0; i < 4; ++i)
		a.step();
	if (s->getApparentRightAscension() > 89. && s->getApparentRightAscension() < 91.)
		std::cout << "\033[1;32mstepApparentRightAscension() loops around as expected with a value of " << s->getApparentRightAscension() << ".\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! stepApparentRightAscension() does not loop around as expected with a value of " << s->getApparentRightAscension() << "!\033[0m" << std::endl;
	
	for (int i = 0; i < 360; ++i)
		a.step();
	if (s->getApparentDeclination() > 23)
		std::cout << "\033[1;32mstepApparentDeclination() works as expected with a value of " << s->getApparentDeclination() << ".\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! stepApparentDeclination() does not perform as expected with a value of " << s->getApparentDeclination() << "!\033[0m" << std::endl;
	
	for (int i = 0; i < 740; ++i)
		a.step();
	if (s->getApparentDeclination() < -23)
		std::cout << "\033[1;32mstepApparentDeclination loops around as expected with a value of " << s->getApparentDeclination() << ".\033[0m" << std::endl;
	else
		std::cout << "\033[1;31mError! stepApparentDeclination does not loop around as expected with a value of " << s->getApparentDeclination() << "!\033[0m" << std::endl;
}

int main(int argc, char** argv){
	std::cout << "\nTesting class \033[1;34m'sunEngine'\033[0m:" << std::endl;
	
	std::cout << "\n\033[1;37mTesting constructors:\033[0m" << std::endl;
	testConstructors();
	
	std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
	testMiscellany();
	
	std::cout << "\nTesting of class \033[1;34m'sunEngine'\033[0m complete." << std::endl;
}