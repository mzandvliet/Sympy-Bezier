import math
import random;
from sympy import *

'''
=== Toying with Benford's law ===

https://en.wikipedia.org/wiki/Benford's_law

Introductory video by Numberphile: www.youtube.com/watch?v=XXjlR2OK1kM
Some example datasets where Benford's law holds: https://testingbenfordslaw.com


'''

def benford_estimate():
    # Benford's function for estimating the prevalence of leading digit

    base = 10 # base 10, decimal

    a = symbols('a')
    benf = log((a+1)/a, base)
    
    for i in range(base):
        print("%s = %s"%(i, benf.evalf(subs={a: 1})))

def benford2():
    numStudents = 9999
    numRolls = 30

    maxDigits = math.ceil(math.log(6 ** numRolls, 10))
    print("max digits: %s"%(maxDigits))
    
    numOnes = 0
    for student in range(numStudents):
        sum = 1
        for roll in range(numRolls):
            value = random.randint(1, 6)
            sum *= value

        print(sum)

        if (is_first_digit_one(sum, maxDigits)):
            print("yes")
            numOnes += 1
    
    print(numOnes / numStudents * 100.0)


def benford3():
    '''

    Imagine we have 60 students, each with pen, paper and dice.

    Each student rolls a 6-sided die a given number of times.
    They multiply each of their rolls together, for example:

    product of dice rolls = 1 * 5 * 4 * 3 * 4 * 3 * 6  * ... * 4 = 472

    When they are done we have a list of 60 large numbers

    The leading digit of any number can be any of [0 ... 9], and so we
    count how many times 1 is the leading number, 2 is the leading number,
    and so on. For example, in the list:

    [
        472,
        102,
        11,
        1,
        30001,
        170
    ]

    we have 4 numbers leading with 1
    we have 1 number  leading with 3
    we have 1 number  leading with 4

    So in this list 4/6 * 100 = 60% of the numbers have 1 as the leading digit

    When I run this program I get:

    0 : 0.0 percent
    1 : 23.333333333333332 percent
    2 : 20.0 percent
    3 : 13.333333333333334 percent
    4 : 16.666666666666664 percent
    5 : 8.333333333333332 percent
    6 : 8.333333333333332 percent
    7 : 5.0 percent
    8 : 1.6666666666666667 percent
    9 : 3.3333333333333335 percent

    When I run this program with a larger amount of students and dice rolls, I get:

    0 : 0.0 percent
    1 : 33.33333333333333 percent
    2 : 18.333333333333332 percent
    3 : 11.666666666666666 percent
    4 : 10.0 percent
    5 : 11.666666666666666 percent
    6 : 10.0 percent
    7 : 1.6666666666666667 percent
    8 : 1.6666666666666667 percent
    9 : 1.6666666666666667 percent

    '''

    numStudents = 60
    numRolls = 15

    # Calculate the maximum amount of digits needed to write
    # the products we'll get.
    maxDigits = math.ceil(math.log(6 ** numRolls, 10))
    print("max digits: %s" % (maxDigits))

    # We will track for each possible leading digits how many
    # times we see it, and store it in this array
    digitCounts = [0] * 10

    # Each student...
    for student in range(numStudents):
        # Rolls their dice and multiplies the results of the rolls together
        product = 1
        for roll in range(numRolls):
            value = random.randint(1, 6)
            product *= value

        # And we keep track of wich value the leading digit of the product had
        digitCounts[first_digit_value(product, maxDigits)] += 1

    # Print the results in the console
    for i in range(0, 10):
        print("%s : %s percent"%(i, digitCounts[i] / numStudents * 100.0))


def benford4():
    '''

    '''

    numStudents = 9999
    numRolls = 30

    maxDigits = math.ceil(math.log(6 * numRolls, 10))
    print("max digits: %s" % (maxDigits))

    digitCounts = [0] * 10
    for student in range(numStudents):
        sum = 0
        for roll in range(numRolls):
            value = random.randint(1, 6)
            sum += value

        # print(sum)

        digitCounts[first_digit_value(sum, maxDigits)] += 1

    for i in range(0, 10):
        print("%s : %s percent" % (i, digitCounts[i] / numStudents * 100.0))


def is_first_digit_one(number, maxDigits):
    '''
    Returns True if the leading digit of number is 1, False if not.

    Parameters:
        number (int): A decimal integer
        maxDigits (int): A natural number indicating the maximum amount of digits to check
    '''

    # Run through all decimal digits from lowest to highest, keep track of which
    # last non-zero number was seen
    result = False
    for d in range(maxDigits):
        digit = get_digit(number, d)
        if (digit == 1):
            result = True
        elif (digit > 1):
            result = False

    return result


def first_digit_value(number, maxDigits):
    '''
    Returns the value of the leading non-zero digit in number

    Parameters:
        number (int): A decimal integer
        maxDigits (int): A natural number indicating the maximum amount of digits to check
    '''

    digit = 0
    for d in range(maxDigits):
        digitValue = get_digit(number, d)
        if (digitValue > 0):
            digit = digitValue

    return digit


def first_digit(number, maxDigits):
    '''
    Returns the index (position) of the leading non-zero digit in number

    Parameters:
        number (int): A decimal integer
        maxDigits (int): A natural number indicating the maximum amount of digits to check
    '''

    digit = 0
    for d in range(maxDigits):
        digitValue = get_digit(number, d)
        if (digitValue > 0):
            digit = d

    return digit

def get_digit(number, n):
    '''
    Returns the value of the n'th digit in number

    Parameters:
        number (int): A decimal integer
        n (int): A natural number indicating the index of the digit in number to check
    '''
    return number // 10**n % 10

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)

    benford3()


if __name__ == "__main__":
    main()
