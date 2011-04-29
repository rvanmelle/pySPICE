
class affnum:
	
	ERROR_TERMS = {}
	
	def __init__(self, x, coef):
		self.X = x
		self.coef = coef
		
	def expand(self):
		pairs = []
		for x in self.coef:
			pairs.append((self.X))
			for y in self.coef
		
	def __str__(self):
		t = "%s" % self.X
		for x in self.coef:
			t += " + %s*%s" % (self.coef[x], x)
		return t
		
	def __add__(self, other):
		new_coef = {}
		new_val = self.X + other.X
		for x in set(self.coef).union(other.coef):
			sum = 0.0
			if self.coef.has_key(x):
				sum += self.coef[x]
			if other.coef.has_key(x):
				sum += other.coef[x]
			new_coef[x] = sum 
				
		return affnum(new_val, new_coef)
			
		
	def __sub__(self, other):
		pass
		
	def __mul__(self, other):
		pass
		
	def __div__(self, other):
		pass
		
if __name__ == "__main__":
	x = affnum(10, {'e3':2.0, 'e8':-6})
	y = affnum(20, {'e4':3.0, 'e8':4})
	z = x + y
	
	print set([1,2,3])
	print x
	print y
	print z