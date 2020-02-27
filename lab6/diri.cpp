int dirichlet(int *ja, int *ia, double *a, double *b, FILE * macierz, FILE * wektor) 
	{
//numeruje niezerowe elementy A
	int k = -1;
//liczba niezerowych el
	int nz_num = 0;

	for(int l = 0; l < N(nx, ny); ++l) 
	{
		int brzeg = 0;  //0 srodek, 1 brzeg
		double vb = 0;  //potencjal na brzegu
		if(i(nx, l) == 0){
			brzeg = 1;
			vb = V1;
		}

		else if(i(nx, l) == nx)
		{
			brzeg = 1;
			vb = V3;
		}

		else if(j(nx, l) == ny)
		{
			brzeg = 1;
			vb = V2;
		}

		else if(j(nx, l) == 0)
		{
			brzeg = 1;
			vb = V4;
		}

        b[l] = (-1)*(ro1(delta*i(nx,l), delta*j(nx,l),xmax,ymax,xmax/10) + ro2(delta*i(nx,l), delta*j(nx,l),xmax,ymax,xmax/10)); //sigma


		if(brzeg == 1)
			b[l] = vb;


		ia[l] = -1;

		if(l - nx - 1 > 0 && brzeg == 0)
		{
			k++;
			if(ia[l] < 0)
                ia[l] = k;

			a[k] = el(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l - nx - 1;
		}
//poddiag
		if(l-1 > 0 && brzeg == 0) 
		{
			k++;
			if(ia[l] < 0)
                ia[l] = k;
			a[k] = el(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l-1;
		}
//diag
		k++;
		if(ia[l] < 0)
            ia[l] = k;

		if(brzeg == 0)
			a[k] = -(2*el(nx,l,eps1,eps2)+el(nx,l+1,eps1,eps2)+el(nx,l+nx+1,eps1,eps2))/(delta*delta);
		else
			a[k] = 1;

		ja[k] = l;
//naddiag
		if(l < N(nx, ny) && brzeg == 0)
		{
			k++;
			a[k] = el(nx,l+1,eps1,eps2)/(delta*delta);
			ja[k] = l + 1;
		}
//prawa skrajna przek
		if(l < N(nx, ny)-nx-1 && brzeg == 0)
		{
			k++;
			a[k] = el(nx,l+nx+1,eps1,eps2)/(delta*delta);
			ja[k] = l + nx + 1;
		}
        if(l%5 == 0 && l != 0)
            fprintf(wektor, "\n");
        fprintf(wektor,"%d %d %d %f \n", l, i(nx,l), j(nx,l), b[l]);

    }

	nz_num = k+1;
	ia[N(nx, ny)] = nz_num;
    for(int z = 0; z < 5*N(nx,ny); z++)
        fprintf(macierz,"%d %0.f \n", z, a[z]);

    return nz_num;
}
