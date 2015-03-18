#include <stdio.h>
#include <matio.h>

int main() {
    printf("%d.%d.%d",MATIO_MAJOR_VERSION,MATIO_MINOR_VERSION,MATIO_RELEASE_LEVEL);
    exit(0);
}
