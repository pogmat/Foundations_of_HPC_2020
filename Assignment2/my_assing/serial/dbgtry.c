#include <stdio.h>

int main(int argc, char** argv)
{
	FILE* fp = fopen(argv[1], "rb");

	unsigned char begin[20];
	unsigned char word[2];
	fread(begin, 1, 19, fp);
	size_t i = 19;
	begin[19] = 0;
	printf("%s\n", begin);
	while(fread(word, 1, 2, fp)) {
		++i;
		printf("%10ld: %x  %c%c\n", i,*word, word[0], word[1]);
	}

	fclose(fp);
	
	return 0;
}
