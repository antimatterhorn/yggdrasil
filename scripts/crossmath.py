''' 
op = A + 13 * B / C + D + 12 * E - F - 11 + G * H / I - 10 = 66

'''

def pemdas(a,b,c,d,e,f,g,h,i):
    return a + (13 * b / c) + d + (12 * e) - f - 11 + (g * h / i) - 10
def orderOp(a,b,c,d,e,f,g,h,i):
    sum = a + 13
    sum = sum * b
    sum = sum / c
    sum = sum + d + 12
    sum = sum * e
    sum = sum - f - 11 + g
    sum = sum * h
    sum = sum / i - 10
    return sum


print("Starting permutations...")

from itertools import permutations
print(permutations)

for perm in permutations(range(1, 10), 9):  # range(1, 10) includes 1–9
    A, B, C, D, E, F, G, H, I = perm
    if pemdas(A, B, C, D, E, F, G, H, I) == 66:
        print(f"pem A={A}, B={B}, C={C}, D={D}, E={E}, F={F}, G={G}, H={H}, I={I} = {pemdas(A, B, C, D, E, F, G, H, I)}")
    if orderOp(A, B, C, D, E, F, G, H, I) == 66:
        print(f"order A={A}, B={B}, C={C}, D={D}, E={E}, F={F}, G={G}, H={H}, I={I} = {orderOp(A, B, C, D, E, F, G, H, I)}")