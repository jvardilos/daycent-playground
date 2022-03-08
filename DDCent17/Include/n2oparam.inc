!****************************************************************************
!**
!**  FILE:    N20param.h
!**
!**    Make parameters for the N2O emission factors addressable and
!**       available for input
!**
!**  AUTHOR:  Kendrick Killian
!**           Natural Resource Ecology Laboratory
!**           Colorado State University
!**           1/2015
!**
!****************************************************************************


!
!  standard nitrification and denitrification parameters
!  cNnh4(2)   = {1.0f, -0.0105f, 0.0f, 0.0f};      - ammonium concentration parameters              fNnh4 = 1.0f - f_exponential(nh4_conc, nfc->cNnh4);           nitrify  unused
!  cNph(4)    = {5.0f, 0.56f, 1.0f, 0.45f};        - ph parameters                                  fNph = f_arctangent(layers->lyrpH[1]*(*pHscale), nfc->cNph);  nitrify
!  cNsoilt(4) = {35.0f, -5.0f, 4.5f, 7.0f};        - temperature effect parameters                  abiotic = max(fNwfps * fNsoilt, sitepar->Ncoeff);             nitrify
!  cNwfps(2)  = {30.0f, -9.0f};                    - Water Filled Pore Space parameters exponential abiotic = max(fNwfps * fNsoilt, sitepar->Ncoeff);             nitrify
!  cDno3(4)   = {9.23f, 1.556f, 76.91f, 0.00222f}; - nitrate parameters from Del Grosso

!typedef struct {
!  float cNph(4);
!  float cNsoilt(4);
!  float cNwfps(2);
!  float cDno3(4);
!  float Nresolim ;
!} NFCON, *NFCON_PT;

   MODULE n2oparam
         USE ISO_C_BINDING
     real             (c_float), save, bind(C, name="cNph") ::   cNph(4);
     real             (c_float), save, bind(C, name="cNsoilt") ::   cNsoilt(4);
     real             (c_float), save, bind(C, name="cNwfps") ::   cNwfps(2);
     real             (c_float), save, bind(C, name="cDno3") ::   cDno3(4);
     real             (c_float), save, bind(C, name="Nresolim") ::   Nresolim;
   END MODULE n2oparam
