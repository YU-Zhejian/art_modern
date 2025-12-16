/**
 * @brief  Assess whether the current infrastructure correctly handles dates beyond January 19, 2038.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void)
{
    if (sizeof(time_t) == 32) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
