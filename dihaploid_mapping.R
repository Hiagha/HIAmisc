r <- 0.3
bivalent <- data.frame(one = c(1,2,3,4), two = c(1,3,2,4), three = c(1,4,2,3))
quadrivalent <- matrix(c(1,2,3,4))
transitional.prob.quad<- matrix(c(1 - r, r/3, r/3, r/3, r/3, 1 - r, r/3, r/3, r/3, r/3, 1 - r, r/3, r/3, r/3, r/3, 1 - r), nrow = 4)
transitional.prob.bi <- matrix(c(1 - r, r, r, 1 - r), nrow = 2)
bivalent.transition <- kronecker(transitional.prob.bi, transitional.prob.bi)
quad.transition <- kronecker(transitional.prob.quad, transitional.prob.quad)
bivalent.by.quad <- kronecker(transitional.prob.quad, transitional.prob.bi)

