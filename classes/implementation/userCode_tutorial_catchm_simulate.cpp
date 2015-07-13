const double f= 2.778e-07; // Factor to convert mm/h into m/s
double a, b;               // Aggregated parameters of the model
double T_new;              // Temporary, to store the new water content

// Compute new soil water content
a= -1. * f * paramNum(Dsat) / paramNum(M) / (paramNum(Tsat) - paramNum(Tmin));
b= inputExt(P) * f / paramNum(M) - a * paramNum(Tmin);
T_new= min( paramNum(Tsat), (stateScal(T) + b/a) * exp(a * delta_t) - b/a);

// Set outputs
//   (1) Average runoff rate computed from water balance
set_output(R)= inputExt(P) * f * paramNum(A) -
               (T_new - stateScal(T)) / delta_t * paramNum(M) * paramNum(A);
//   (2) Soil water content at end of time step
set_output(Tout)= T_new;

// Update soil water content
set_stateScal(T)= T_new;

