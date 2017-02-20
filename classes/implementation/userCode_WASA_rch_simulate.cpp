
double qin = inputSim(q_in);
double qout = -9999.;

//////////////////////////////////////////////////////////
// ROUTING
// Calculates maximum river reach outflow (m3/s)
//////////////////////////////////////////////////////////
	
// initialise unit hydrograph vector object from paramFun input
vector<double> uh_f(paramNum(nuh));
for (unsigned int i=0; i<uh_f.size(); i++) {
	uh_f[i] = paramFun(uh, i+1);
}
// build object of type vector to pass uh_q to function chan_route()
vector<double> uh_remain(paramNum(nuh));
uh_remain = stateVect(uh_q);

// calculate river reach output from external function
qout = chan_route(qin, uh_f, uh_remain, sharedParamNum(choice_route));



//////////////////////////////////////////////////////////
// TRANSMISSION LOSSES (m3/s)
// Substracted from river output.
//////////////////////////////////////////////////////////

if (qout > 1e-6) {
	double loss = trans_loss(qout, inputSim(pet), paramNum(chan_len), sharedParamNum(choice_transloss));
	qout -= loss;
}

//////////////////////////////////////////////////////////
// OUTPUT
// Actual river reach outflow (m3/s)
//////////////////////////////////////////////////////////

set_output(q_out) = qout;

// update state vector
set_stateVect(uh_q) = uh_remain;
