def add_assumptions1(k):
  if k==0:
    return
  assume(x<=k)
  print k, x-k, assumptions()
  add_assumptions1(k-1)
  #print "ahoj"
  print k, x-k, assumptions() 
  forget(x<=k)

def add_assumptions2(k):
  if k==0:
    return 
  assume(x-k<=0)
  print k, x-k, assumptions() 
  add_assumptions2(k-1)
  print k, x-k, assumptions() 
  forget(x-k<=0)

  
# sage: add_assumptions1(3)
# 3 x - 3 [x <= 3]
# 2 x - 2 [x <= 3, x <= 2]
# 1 x - 1 [x <= 3, x <= 2, x <= 1]
# sage: add_assumptions2(3)
# 3 x - 3 [x - 3 <= 0]
# 2 x - 2 [x - 3 <= 0, x - 2 <= 0]
# 1 x - 1 [x - 3 <= 0, x - 2 <= 0]
