#include <stdio.h>
int main() {
	char stringa[5] = "prova";
	int i = 1;
	switch (i){
		case 1:
			printf("1 %d\n",i);
		case 2:
			printf("2 %d\n",i);
		case 3:
			printf("3 %d\n",i);
		default:
			printf("altro\n");
	}

	if (i==3) {
		printf("3 %d\n",i);
	}
	if (i==2) {
		printf("2 %d\n",i);
	} else if (i==1) {
		printf("1 %d\n",i);
	} else {
		printf("altro\n");
	}



	return 0;
}