import random


if __name__ == "__main__":
    rng = random.SystemRandom()
    print("".join(rng.choices("ACGT", k=150)))
