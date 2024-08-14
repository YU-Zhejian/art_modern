__inline__ static int f(int i)
{
    return i - 1;
}

int main(void)
{
    volatile int j = f(1);
    return j;
}
