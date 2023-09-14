% x*log(x/y)

function z = xlogxdy(x,y) 

z = x*log((x+eps)/(y+eps)); 