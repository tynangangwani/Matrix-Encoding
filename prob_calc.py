import math
def calc_prob(n,k, t):
   return ((1/65537)**(n-k))*sum([math.comb(n,i)*((65536)**i) for i in range(t)])
print(calc_prob(16,8,5))
print(calc_prob(16,8,5)**12)

print(calc_prob(150,140,6))
print(calc_prob(150,140,6)**12)
