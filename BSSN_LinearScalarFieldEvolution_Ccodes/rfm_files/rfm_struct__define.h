for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct.f0_of_xx0[i0] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1.0/SINHW) - exp(-1/SINHW));
}

for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct.f0_of_xx0__D0[i0] = AMPL*(exp(xx0/SINHW)/SINHW + exp(-xx0/SINHW)/SINHW)/(exp(1.0/SINHW) - exp(-1/SINHW));
}

for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct.f0_of_xx0__DD00[i0] = AMPL*(exp(xx0/SINHW)/pow(SINHW, 2) - exp(-xx0/SINHW)/pow(SINHW, 2))/(exp(1.0/SINHW) - exp(-1/SINHW));
}

for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct.f0_of_xx0__DDD000[i0] = AMPL*(exp(xx0/SINHW)/pow(SINHW, 3) + exp(-xx0/SINHW)/pow(SINHW, 3))/(exp(1.0/SINHW) - exp(-1/SINHW));
}

for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
    const REAL xx1 = xx[1][i1];
    rfmstruct.f1_of_xx1[i1] = sin(xx1);
}

for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
    const REAL xx1 = xx[1][i1];
    rfmstruct.f1_of_xx1__D1[i1] = cos(xx1);
}

for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
    const REAL xx1 = xx[1][i1];
    rfmstruct.f1_of_xx1__DD11[i1] = -sin(xx1);
}

