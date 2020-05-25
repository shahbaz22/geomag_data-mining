a=True
b=False
def function(c):
	if c:
		print('its true')
	else:
		print('its false')
function(a)
# The Boolean expression after the if statement is called the condition. 
# If it is true, then all the indented statements get executed. 
# If not, then all the statements indented under the else clause get executed.
# if BOOLEAN EXPRESSION:
#     STATEMENTS_1        # Executed if condition evaluates to True
# else:
#     STATEMENTS_2        # Executed if condition evaluates to False
if True:          # This is always True,
    pass          #   so this is always executed, but it does nothing
else:
    pass

name = input("Enter your name")
name= 'Hello ' + name
print(name)
