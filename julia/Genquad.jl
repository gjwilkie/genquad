module Genquad
using Cubature
using Base.Test

export genquadrules

function genquadrules(wgtfcn::Function,N::Integer,endpts::Array{Float64,1})
   functest = wgtfcn(0.5*(endpts[1]+endpts[2]))
   intwgts = zeros(typeof(functest),N)
   abscissae = zeros(typeof(endpts[1]),N)
   const epsfac = 10.0

   a=zeros(typeof(functest),N)
   b=zeros(typeof(functest),N)

   epsilon = eps(typeof(functest))
   order = zero(typeof(N))

   # Sanity check on endpts
   length(endpts) >= 2            || error("Must pass at least two end points to genquadrules")
   for i in 1:length(endpts)-1
      endpts[i] <  endpts[i+1]   || error("endpts must be monotonically increasing")
   end


   function wgt_x(x::Real)
      return x*wgtfcn(x)
   end

   function poly_sq(x::Real)  
      ym1 = zero(typeof(x))
      y= one(typeof(x))
      for i = 1:order
         temp = y 
         y = (x - a[i])*y - b[i]*ym1
         ym1 = temp
      end
      return wgtfcn(x)*y*y
   end

   function poly_sqx(x::Real)  
      ym1 = zero(typeof(x))
      y= one(typeof(x))
      for i = 1:order
         temp = y 
         y = (x - a[i])*y - b[i]*ym1
         ym1 = temp
      end
      return wgtfcn(x)*y*y*x
   end

   moment1, dummy = pquad(wgt_x,endpts...;reltol=epsfac*epsilon)
   moment2, dummy = pquad(wgtfcn,endpts...;reltol=epsfac*epsilon)
   b[1] = moment2
#   b[1] = 1.0
   a[1] = moment1/moment2
   order = 1
   moment1, dummy = pquad(poly_sqx,endpts...;reltol=epsfac*epsilon)
   moment2, dummy = pquad(poly_sq,endpts...;reltol=epsfac*epsilon)
   a[2] = moment1/moment2
   order = 0
   moment1, dummy = pquad(poly_sq,endpts...;reltol=epsfac*epsilon)
   b[2] = moment2/moment1
   for n = 3:N
      order = n-1
      moment1, dummy = pquad(poly_sqx,endpts...;reltol=epsfac*epsilon)
      moment2, dummy = pquad(poly_sq,endpts...;reltol=epsfac*epsilon)
      a[n] = moment1/moment2
      order = n-2
      moment1, dummy = pquad(poly_sq,endpts...;reltol=epsfac*epsilon)
      b[n] = moment2/moment1
   end

   b = copy(sqrt(b))

#println("a=")
#println(a)
#println("b=")
#println(b)

   matrix = zeros(N,N)
   for i = 1:N
      matrix[i,i] = a[i]
      if i > 1
         matrix[i,i-1] = b[i] 
      end
      if i < N
         matrix[i,i+1] = b[i+1] 
      end
   end

   vals,vecs = eig(matrix)
   
   for i = 1:N
      abscissae[i] = vals[i]
      intwgts[i] = (b[1]*vecs[1,i])^2
   end

   return abscissae, intwgts
end

"Acts as a wrapper for Cubature's pquadrature, with the possibility of infinite or semiinfinite domains"
function pquad(f::Function,endpts...;reltol=1.0e-8)

#   f_y = Function
   endpts_y = [0.0,0.0]
   # If user specifies any singularities in the domain, just use quadgk
   if length(endpts) > 2
      return Base.quadgk(f,endpts...;reltol=reltol)
   else
      if isinf(endpts[1])
         if isinf(endpts[2])
            # If infinite:  
            f_y(x) = f(x/(1.0-x^2))*(1+x^2)/(1.0-x^2)^2 
            endpts_y = [-1.0,1.0]
         else  
            # If semi-infinite:  
            f_y(x) = -f(endpts[2]-x/(1.0-x))/(1.0-x)^2
            endpts_y = [0.0,1.0]
         end
      else
         if isinf(endpts[2])
#println("Got here")
            # If semi-infinite:  
            f_y(x) = f(endpts[1] + x/(1.0-x))/(1.0-x)^2
            endpts_y = [0.0,1.0]
         else  
            # If both bounds are finite, do no transformation
            f_y(x) = f(x)
            endpts_y = endpts
         end
      end
      if (isinf(f_y(endpts_y[1])) || isinf(f_y(endpts_y[2])))  || (isnan(f_y(endpts_y[1])) || isnan(f_y(endpts_y[2]))) 
         # If endpoints cannot be evaluated, just use quadgk
         return Base.quadgk(f,endpts...;reltol=reltol)
      else
         return Cubature.hquadrature(f_y,endpts_y[1],endpts_y[2];reltol=reltol)
      end
   end

   # Transformed integrand
   function fy(x) 
      return f(a*x+b)
   end
end

function test()

   const epsfac = 10.0
   Nmax = 60
   tol = 100.0*eps(Float64)*epsfac

   function testpoly(x)
      return x.^9 - 2.0*x.^8 + 11.0*x.^7 + 50.0*x.^6 - 60.0*x.^5 - 40.0*x.^4 + 19.0*x.^3 - 7.0*x.^2 + 2.0*x - 5.0
   end

   function testtrig(x)
      return cos(4.0*x - 2.0)
#      return tanh(4.0*x - 2.0)
   end

   # Hermite quadrature
   print("Testing Hermite quadrature...")
   function g(x)
      return exp(-x.^2)
   end

   function gtrig(x) 
      return testtrig(x)*g(x)
   end
   function gpoly(x) 
      return testpoly(x)*g(x)
   end
   
   xpts,wgts = genquadrules(g,5,[-Inf,Inf])
   ans,dummy = quadgk(gpoly,-Inf,Inf,reltol=eps(Float64)*epsfac)
   guess = sum(wgts.*testpoly(xpts))
   abs(ans - guess)/abs(ans) <= tol || error("N=5 should exactly recover polynomial of degree 9, but is off by "string(ans-guess)" with an exact answer of "string(ans)" +/- "string(dummy)". Tol = "string(tol))
   
   converged = false
   n=4
   while !converged
      n <= Nmax || error("Could not converge integrals within "string(Nmax)" test points.")
      
      xpts,wgts = genquadrules(g,n,[-Inf,Inf])
      ans,dummy = quadgk(gtrig,-Inf,Inf,reltol=eps(Float64)*epsfac)
      guess = sum(wgts.*testtrig(xpts))
      if (guess - ans) <= tol
         converged = true
      end
      n += 1
   end
   println(" PASSED. Converged to exact answer at n = "string(n-1))

   # Shifted Laguerre quadrature
   print("Testing shifted Laguerre quadrature...")
   function g(x)
      return exp(-x)
   end

   function gtrig(x) 
      return testtrig(x)*g(x)
   end
   function gpoly(x) 
      return testpoly(x)*g(x)
   end
   
   xpts,wgts = genquadrules(g,5,[-4.0,Inf])
   ans,dummy = quadgk(gpoly,-4.0,Inf,reltol=eps(Float64)*epsfac)
   guess = sum(wgts.*testpoly(xpts))
   abs(ans - guess)/abs(ans) <= tol || error("N=5 should exactly recover polynomial of degree 9, but is off by "string(ans-guess)" with an exact answer of "string(ans)" +/- "string(dummy))
   
   converged = false
   n=4
   while !converged
      n <= Nmax || error("Could not converge integrals within "string(Nmax)" test points.")
      
      xpts,wgts = genquadrules(g,n,[-4.0,Inf])
      ans,dummy = quadgk(gtrig,-4.0,Inf)
      guess = sum(wgts.*testtrig(xpts))
      if (guess - ans) <= tol
         converged = true
      end
      n += 1
   end
   println(" PASSED. Converged to exact answer at n = "string(n-1))
 
   # Shifted Legendre quadrature
   print("Testing shifted Legendre quadrature...")
   function g(x)
      return 4.0
   end

   function gtrig(x) 
      return testtrig(x)*g(x)
   end
   function gpoly(x) 
      return testpoly(x)*g(x)
   end
   
   xpts,wgts = genquadrules(g,5,[3.0,20.0])
   ans,dummy = quadgk(gpoly,3.0,20.0,reltol=eps(Float64)*epsfac)
   guess = sum(wgts.*testpoly(xpts))
   abs(ans - guess)/abs(ans) <= tol || error("N=5 should exactly recover polynomial of degree 9, but is off by "string(ans-guess)" with an exact answer of "string(ans)" +/- "string(dummy))
   
   converged = false
   n=4
   while !converged
      n <= Nmax || error("Could not converge integrals within "string(Nmax)" test points.")
      
      xpts,wgts = genquadrules(g,n,[3.0,20.0])
      ans,dummy = quadgk(gtrig,3.0,20.0)
      guess = sum(wgts.*testtrig(xpts))
      if (guess - ans) <= tol
         converged = true
      end
      n += 1
   end
   println(" PASSED. Converged to exact answer at n = "string(n-1))
 
end

end
