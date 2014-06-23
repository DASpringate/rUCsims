\name{lin_sim}
\alias{lin_sim}
\title{This function generates the simulated data given the parameters
Return is an object of class 'UCsims' which is a list containing a dataframe of the simulated data and a vector of parameter values.
For generality, the confounders are labelled U1 and U2.
In this case, U1 is binary and U2 is continuous}
\usage{
lin_sim(obs = 1000, orUB = 0.1, probE_0 = 0.2, probY_0 = 0.7,
  orEY_0 = 1, probU = 0.4, modUB = 1, orME = 2, orMY_0 = 2, rho = 0,
  orUBe = NULL, orUBy_0 = NULL)
}
\description{
This function generates the simulated data given the
parameters Return is an object of class 'UCsims' which is a
list containing a dataframe of the simulated data and a
vector of parameter values. For generality, the confounders
are labelled U1 and U2. In this case, U1 is binary and U2
is continuous
}

