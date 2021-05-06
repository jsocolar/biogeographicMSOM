parameters{
  real mu;
}
model{
  mu ~ logistic(0,1);
}
generated quantities{
  real p = inv_logit(mu);
}