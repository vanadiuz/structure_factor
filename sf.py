def calculateStructureFactor(x, y, z, L, order, sepz=False, res_xyz=1, res_z=1):
   ''' ACHTUNG! res_xyz also for xy-plane'''
   from math import pi,sin,cos,sqrt
   order2=order*order
   twoPI_L=2*pi/L
   mult_coeff = res_xyz
   ff=[0.0 for i in range(2*order2*mult_coeff**2)]
   if not sepz:
      for i_new in range((order+1)*mult_coeff):
         i = i_new/mult_coeff 
         for l_new in range((2*order+1)*mult_coeff):
            j=-order+l_new/mult_coeff 
            j_new = -order*mult_coeff+l_new
            for m_new in range((2*order+1)*mult_coeff):
               k=-order+m_new/mult_coeff
               k_new=-order*mult_coeff+m_new
               n=i*i+j*j+k*k
               n_arr = i_new*i_new + j_new*j_new + k_new*k_new
               if ((n<=order2) & (n>=1)):
                  C_sum=0.
                  S_sum=0.
                  for p in range(len(x)):
                     qr = twoPI_L * (i*x[p] + j*y[p] + k*z[p])
                     C_sum=C_sum+cos(qr)
                     S_sum=S_sum+sin(qr)
                  ff[2*n_arr-2]=ff[2*n_arr-2]+C_sum*C_sum + S_sum*S_sum
                  ff[2*n_arr-1]=ff[2*n_arr-1]+1
      n=len(x)
      for i in range(order2*mult_coeff**2):
         if ff[2*i+1]!=0:
            ff[2*i]=ff[2*i]/(n*ff[2*i+1])

      qs=[]
      Sq=[]
      for i in range(order2*mult_coeff**2):
         if ff[2*i+1]>0.:
            qs.append(twoPI_L*sqrt(i+1.)/mult_coeff)
            Sq.append(ff[2*i])
      return qs,Sq
   else:
      k=0
      for i_new in range((order+1)*mult_coeff):
         i = i_new/mult_coeff 
         for l_new in range((2*order+1)*mult_coeff):
            j=-order+l_new/mult_coeff 
            j_new = -order*mult_coeff+l_new
            n=i*i+j*j
            n_arr = i_new*i_new + j_new*j_new
            if ((n<=order2) & (n>=1)):
               C_sum=0.
               S_sum=0.
               for p in range(len(x)):
                  qr = twoPI_L * (i*x[p] + j*z[p])
                  C_sum=C_sum+cos(qr)
                  S_sum=S_sum+sin(qr)
               ff[2*n_arr-2]=ff[2*n_arr-2]+C_sum*C_sum + S_sum*S_sum
               ff[2*n_arr-1]=ff[2*n_arr-1]+1
      n=len(x)
      for i in range(order2*mult_coeff**2):
         if ff[2*i+1]!=0:
            ff[2*i]=ff[2*i]/(n*ff[2*i+1])

      qsxy=[]
      Sqxy=[]
      for i in range(order2*mult_coeff**2):
         if ff[2*i+1]>0.:
            qsxy.append(twoPI_L*sqrt(i+1.)/mult_coeff)
            Sqxy.append(ff[2*i])

      mult_coeff = res_z
      ff=[0.0 for i in range(2*order2*mult_coeff**2)]
      for k_new in range((order+1)*mult_coeff):
         k = k_new/mult_coeff
         n=k*k
         n_arr = k_new*k_new
         if ((n<=order2) & (n>=1)):
            C_sum=0.
            S_sum=0.
            for p in range(len(y)):
               qr = twoPI_L * (k*y[p])
               C_sum=C_sum+cos(qr)
               S_sum=S_sum+sin(qr)
            ff[2*n_arr-2]=ff[2*n_arr-2]+C_sum*C_sum + S_sum*S_sum
            ff[2*n_arr-1]=ff[2*n_arr-1]+1
      n=len(y)
      for i in range(order2*mult_coeff**2):
         if ff[2*i+1]!=0:
            ff[2*i]=ff[2*i]/(n*ff[2*i+1])

      qsz=[]
      Sqz=[]
      for i in range(order2*mult_coeff**2):
         if ff[2*i+1]>0.:
            qsz.append(twoPI_L*sqrt((i+1)/mult_coeff**2))
            Sqz.append(ff[2*i])
      return qsxy, Sqxy, qsz, Sqz