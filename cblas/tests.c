  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -6699 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 1)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -81.2 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 2)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -6699 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 3)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -81.2 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 4)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -6699 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 5)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -81.2 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 6)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -6699 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 7)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 82.5 };
   float X[] = { -81.2 };
   int incX = 1;
   float x_expected[] = { -81.2 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 8)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -5495 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 9)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -70 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 10)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -5495 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 11)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -70 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 12)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -5495 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 13)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -70 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 14)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -5495 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 15)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 78.5 };
   double X[] = { -70 };
   int incX = 1;
   double x_expected[] = { -70 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 16)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -194.4, 2913.48 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 17) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 17) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -39.2, -87 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 18) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 18) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -194.4, 2913.48 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 19) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 19) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -39.2, -87 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 20) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 20) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -194.4, 2913.48 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 21) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 21) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -39.2, -87 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 22) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 22) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -194.4, 2913.48 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 23) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 23) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -27, -14.4 };
   float X[] = { -39.2, -87 };
   int incX = 1;
   float x_expected[] = { -39.2, -87 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 24) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 24) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 2987.54, 11951.72 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 25) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 25) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 94.9, -87.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 26) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 26) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 2987.54, 11951.72 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 27) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 27) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 94.9, -87.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 28) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 28) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 2987.54, 11951.72 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 29) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 29) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 94.9, -87.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 30) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 30) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 2987.54, 11951.72 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 31) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 31) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -45.7, 83.9 };
   double X[] = { 94.9, -87.3 };
   int incX = 1;
   double x_expected[] = { 94.9, -87.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 32) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 32) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { 1289.04 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 33)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { -78.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 34)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { 1289.04 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 35)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { -78.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 36)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { 1289.04 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 37)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { -78.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 38)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { 1289.04 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 39)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -16.4 };
   float X[] = { -78.6 };
   int incX = -1;
   float x_expected[] = { -78.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 40)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { 145.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 41)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { -16.5 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 42)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { 145.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 43)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { -16.5 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 44)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { 145.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 45)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { -16.5 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 46)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { 145.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 47)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -8.8 };
   double X[] = { -16.5 };
   int incX = -1;
   double x_expected[] = { -16.5 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 48)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -8727.03, -7094.49 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 49) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 49) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -85.8, -89.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 50) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 50) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -8727.03, -7094.49 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 51) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 51) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -85.8, -89.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 52) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 52) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -8727.03, -7094.49 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 53) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 53) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -85.8, -89.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 54) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 54) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -8727.03, -7094.49 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 55) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 55) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 89.9, -11.3 };
   float X[] = { -85.8, -89.7 };
   int incX = -1;
   float x_expected[] = { -85.8, -89.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 56) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 56) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 1761.99, 3667.07 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 57) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 57) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 77.6, 98.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 58) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 58) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 1761.99, 3667.07 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 59) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 59) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 77.6, 98.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 60) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 60) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 1761.99, 3667.07 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 61) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 61) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 77.6, 98.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 62) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 62) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 1761.99, 3667.07 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 63) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 63) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 31.7, 7.1 };
   double X[] = { 77.6, 98.3 };
   int incX = -1;
   double x_expected[] = { 77.6, 98.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 64) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 64) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { -4398.72 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 65)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { 69.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 66)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { -4398.72 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 67)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { 69.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 68)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { -4398.72 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 69)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { 69.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 70)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { -4398.72 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 71)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -63.2 };
   float X[] = { 69.6 };
   int incX = -1;
   float x_expected[] = { 69.6 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 72)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { -146.16 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 73)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { 23.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 74)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { -146.16 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 75)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { 23.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 76)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { -146.16 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 77)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { 23.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 78)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { -146.16 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 79)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -6.3 };
   double X[] = { 23.2 };
   int incX = -1;
   double x_expected[] = { 23.2 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 80)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -1855.98, -1391.04 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 81) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 81) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -12.6, 21.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 82) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 82) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -1855.98, -1391.04 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 83) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 83) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -12.6, 21.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 84) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 84) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -1855.98, -1391.04 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 85) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 85) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -12.6, 21.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 86) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 86) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -1855.98, -1391.04 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 87) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 87) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -10.8, 91.8 };
   float X[] = { -12.6, 21.7 };
   int incX = -1;
   float x_expected[] = { -12.6, 21.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 88) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 88) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 1770.72, 5926.46 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 89) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 89) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 60.9, -14.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 90) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 90) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 1770.72, 5926.46 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 91) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 91) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 60.9, -14.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 92) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 92) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 1770.72, 5926.46 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 93) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 93) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 60.9, -14.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 94) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 94) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 1770.72, 5926.46 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 95) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 95) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 5.9, 98.7 };
   double X[] = { 60.9, -14.3 };
   int incX = -1;
   double x_expected[] = { 60.9, -14.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 96) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 96) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { 572.73, -5852.49 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 97)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { 906.9, 89.9 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 98)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { -326.27, -5839.06 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 99)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { 7.9, 103.33 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 100)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { -173.44, -5852.49 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 101)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { 160.73, 89.9 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 102)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { -326.27, -5773.49 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 103)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -41.3, 10, 1.7, -65.1 };
   float X[] = { 7.9, 89.9 };
   int incX = 1;
   float x_expected[] = { 7.9, 168.9 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 104)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { -55.72, -826.88 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 105)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { -321.14, -32.3 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 106)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { 260.82, -530.18 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 107)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { -4.6, 264.4 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 108)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { 2344.17, -826.88 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 109)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { 2078.75, -32.3 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 110)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { 260.82, -871.96 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 111)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -56.7, 9.8, -64.5, 25.6 };
   double X[] = { -4.6, -32.3 };
   int incX = 1;
   double x_expected[] = { -4.6, -77.38 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 112)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -13786.78, -614.3, -4074.04, -190.72 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 113) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 113) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -2179.96, 1705.92, 15.2, -38.8 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 114) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 114) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -11662.22, -2404.62, -8192.68, 8381.9 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 115) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 115) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -55.4, -84.4, -4103.44, 8533.82 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 116) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 116) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -15532.1, -1745.58, -4074.04, -190.72 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 117) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 117) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -3925.28, 574.64, 15.2, -38.8 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 118) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 118) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -11662.22, -2404.62, -3511.52, 6516.84 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 119) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 119) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 83.3, -83.5, -58.6, -31.8, -48.6, -80.7, -31.4, -92.7 };
   float X[] = { -55.4, -84.4, 15.2, -38.8 };
   int incX = 1;
   float x_expected[] = { -55.4, -84.4, 577.72, 6668.76 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 120) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 120) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { -337.88, -344.7, -3339.36, 2638.34 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 121) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 121) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { -545.54, -6008.78, 39.2, 56.2 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 122) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 122) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { 239.36, 5585.18, -5570.55, -6587.83 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 123) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 123) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { 31.7, -78.9, -2191.99, -9169.97 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 124) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 124) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { 7444.4, 8153.6, -3339.36, 2638.34 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 125) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 125) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { 7236.74, 2489.52, 39.2, 56.2 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 126) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 126) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { 239.36, 5585.18, -9103.36, 7268.54 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 127) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 127) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -59.9, 27.1, -75.8, -42.6, 90.9, -64.8, 3.7, 62 };
   double X[] = { 31.7, -78.9, 39.2, 56.2 };
   int incX = 1;
   double x_expected[] = { 31.7, -78.9, -5724.8, 4686.4 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 128) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 128) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { -6649.86, 1707.96 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 129)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { 841.14, 17.2 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 130)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { -7573.5, -5106.54 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 131)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { -82.5, -6797.3 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 132)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { -6152.78, 1707.96 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 133)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { 1338.22, 17.2 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 134)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { -7573.5, -2722.29 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 135)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 91.8, 53.7, 82.6, 99.3 };
   float X[] = { -82.5, 17.2 };
   int incX = 1;
   float x_expected[] = { -82.5, -4413.05 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 136)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { 3008.43, -296.07 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 137)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { 24.69, 7.1 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 138)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { 2916.84, -416.49 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 139)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { -66.9, -113.32 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 140)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { 2929.62, -296.07 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 141)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { -54.12, 7.1 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 142)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { 2916.84, -1159.08 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 143)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -43.6, 12.9, 1.8, -41.7 };
   double X[] = { -66.9, 7.1 };
   int incX = 1;
   double x_expected[] = { -66.9, -855.91 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 144)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { -7567.58, -3959.7, 999.08, 6732.94 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 145) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 145) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { -5524.97, -2090.01, 77.1, -39.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 146) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 146) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { -1988.11, -1818.99, 1321.26, 1654.58 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 147) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 147) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { 54.5, 50.7, 399.28, -5118.06 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 148) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 148) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { -7426.67, -4178.57, 999.08, 6732.94 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 149) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 149) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { -5384.06, -2308.88, 77.1, -39.7 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 150) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 150) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { -1988.11, -1818.99, 1103.51, 1604.51 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 151) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 151) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -36.2, 0.3, -45.9, -51.4, -43.3, -52.9, -25.3, 74.3 };
   float X[] = { 54.5, 50.7, 77.1, -39.7 };
   int incX = 1;
   float x_expected[] = { 54.5, 50.7, 181.53, -5168.13 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 152) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 152) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { -5085.89, 4737.18, -53.76, 110.42 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 153) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 153) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { -342.69, 120.78, 2.7, 1.6 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 154) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 154) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { -4702.8, 4681.1, 4675.45, -1979.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 155) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 155) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { 40.4, 64.7, 4731.91, -2088.12 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 156) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 156) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { -4569.52, 4515.29, -53.76, 110.42 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 157) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 157) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { 173.68, -101.11, 2.7, 1.6 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 158) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 158) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { -4702.8, 4681.1, -8948.84, -2959.27 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 159) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 159) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 19.4, 84.8, -95.9, 77.6, 9.6, -67.1, 3.2, 39 };
   double X[] = { 40.4, 64.7, 2.7, 1.6 };
   int incX = 1;
   double x_expected[] = { 40.4, 64.7, -8892.38, -3068.09 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 160) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 160) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { -511.5, 2282.52 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 161)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { -33, -2868.4 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 162)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { -6645.42, 5249.22 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 163)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { -6166.92, 98.3 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 164)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { -511.5, 7308.42 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 165)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { -33, 2157.5 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 166)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { 8325.67, 5249.22 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 167)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 53.4, 89.9, -62.4, 15.5 };
   float X[] = { -33, 98.3 };
   int incX = -1;
   float x_expected[] = { 8804.17, 98.3 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 168)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { 13.98, -6739.9 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 169)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { 0.2, -90.28 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 170)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { -8256.6, -6741.72 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 171)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { -8270.38, -92.1 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 172)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { 13.98, -6723.76 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 173)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { 0.2, -74.14 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 174)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { -824.13, -6741.72 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 175)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 73.2, 9.1, 89.8, 69.9 };
   double X[] = { 0.2, -92.1 };
   int incX = -1;
   double x_expected[] = { -837.91, -92.1 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 176)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { -2849.85, -525.9, 5694.84, 779.02 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 177) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 177) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { 89.8, 22.9, 2357.92, 2729.21 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 178) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 178) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { -4814.55, 973.5, 3374.52, -1935.89 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 179) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 179) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { -1874.9, 1522.3, 37.6, 14.3 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 180) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 180) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { -2849.85, -525.9, -730.38, 2009.71 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 181) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 181) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { 89.8, 22.9, -4067.3, 3959.9 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 182) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 182) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { -1982.91, 759.27, 3374.52, -1935.89 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 183) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 183) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 61.3, -74.8, 31.5, 22.2, -32.4, 52.2, -31.2, 2.1 };
   float X[] = { 89.8, 22.9, 37.6, 14.3 };
   int incX = -1;
   float x_expected[] = { 956.74, 1308.07, 37.6, 14.3 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 184) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 184) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { -8951.38, -4252.82, 4437.45, 10724.69 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 185) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 185) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { 31, 83.1, 4169.11, 7275.79 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 186) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 186) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { -9304.28, -5130.98, 188.64, 3420.7 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 187) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 187) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { -321.9, -795.06, -79.7, -28.2 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 188) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 188) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { -8951.38, -4252.82, -280, 4296.04 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 189) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 189) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { 31, 83.1, -548.34, 847.14 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 190) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 190) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { -16889.23, -5617.63, 188.64, 3420.7 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 191) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 191) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -15.6, -37.4, 93.9, -16.1, 7.4, 8.4, -80.2, 77.8 };
   double X[] = { 31, 83.1, -79.7, -28.2 };
   int incX = -1;
   double x_expected[] = { -7906.85, -1281.71, -79.7, -28.2 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 192) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 192) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { 1671.64, -3972.68 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 193)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { -52.9, -4824.04 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 194)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { 504.15, 820.06 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 195)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { -1220.39, -31.3 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 196)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { 1671.64, -1153.11 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 197)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { -52.9, -2004.47 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 198)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { -1164.14, 820.06 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 199)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { -26.2, 90.6, 37.3, -31.6 };
   float X[] = { -52.9, -31.3 };
   int incX = -1;
   float x_expected[] = { -2888.68, -31.3 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 200)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { -3292.84, -3244.79 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 201)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { 76.4, -3933.34 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 202)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { -4277.32, 659.25 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 203)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { -908.08, -29.3 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 204)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { -3292.84, 3226.29 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 205)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { 76.4, 2537.74 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 206)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { -1795.61, 659.25 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 207)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { -22.5, -51.1, 33.6, -43.1 };
   double X[] = { 76.4, -29.3 };
   int incX = -1;
   double x_expected[] = { 1573.63, -29.3 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 208)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -1784.69, 5174.55, 4626.71, -4971.73 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 209) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 209) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -58.5, 14.6, 1676.95, 2120.49 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 210) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 210) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -9749.6, 3273.42, 3035.06, -7091.32 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 211) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 211) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -8023.41, -1886.53, 85.3, 0.9 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 212) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 212) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -1784.69, 5174.55, 8821.64, -7211.83 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 213) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 213) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -58.5, 14.6, 5871.88, -119.61 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 214) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 214) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -3206.87, 1704.51, 3035.06, -7091.32 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 215) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 215) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   float A[] = { 34.7, -83.5, -17.1, -40.5, -93.6, -21.3, 49.5, -76.1 };
   float X[] = { -58.5, 14.6, 85.3, 0.9 };
   int incX = -1;
   float x_expected[] = { -1480.68, -3455.44, 85.3, 0.9 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 216) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 216) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { 617.92, -4802.36, -4178.81, -6789.39 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 217) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 217) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { 38.9, -42.9, -3675.11, -2672.69 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 218) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 218) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { -87.62, -7065.13, -527.3, -4135 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 219) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 219) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { -666.64, -2305.67, -23.6, -18.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 220) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 220) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { 617.92, -4802.36, 3952.75, -5161.73 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 221) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 221) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { 38.9, -42.9, 4456.45, -1045.03 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 222) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 222) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { -602.09, -2819.64, -527.3, -4135 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 223) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 223) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 2;
   int lda = 2;
   double A[] = { 98.8, 98.6, -8.4, -77.5, 65.1, 45.4, 68.6, -47.8 };
   double X[] = { 38.9, -42.9, -23.6, -18.3 };
   int incX = -1;
   double x_expected[] = { -1181.11, 1939.82, -23.6, -18.3 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 224) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 224) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 3;
   int lda = 3;
   float A[] = { -58, -95.6, 44.3, 61.3, 50.3, 62.9, -41.9, 39.8, -65.6 };
   float X[] = { 73.6, -31.8, 33.1 };
   int incX = 1;
   float x_expected[] = { 237.61, 482.45, -2171.36 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 225)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 3;
   int lda = 3;
   float A[] = { -58, -95.6, 44.3, 61.3, 50.3, 62.9, -41.9, 39.8, -65.6 };
   float X[] = { 73.6, -31.8, 33.1 };
   int incX = 1;
   float x_expected[] = { 4580.01, 2050.19, 33.1 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 226)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 131;
   int N = 3;
   int lda = 3;
   float A[] = { -58, -95.6, 44.3, 61.3, 50.3, 62.9, -41.9, 39.8, -65.6 };
   float X[] = { 73.6, -31.8, 33.1 };
   int incX = 1;
   float x_expected[] = { -4268.8, 2912.14, -6520.84 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 227)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 101;
   int trans = 111;
   int diag = 132;
   int N = 3;
   int lda = 3;
   float A[] = { -58, -95.6, 44.3, 61.3, 50.3, 62.9, -41.9, 39.8, -65.6 };
   float X[] = { 73.6, -31.8, 33.1 };
   int incX = 1;
   float x_expected[] = { 73.6, 4479.88, -4316.38 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 228)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 131;
   int N = 3;
   int lda = 3;
   float A[] = { -58, -95.6, 44.3, 61.3, 50.3, 62.9, -41.9, 39.8, -65.6 };
   float X[] = { 73.6, -31.8, 33.1 };
   int incX = 1;
   float x_expected[] = { -7605.03, -282.16, -2171.36 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 229)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 102;
   int trans = 111;
   int diag = 132;
   int N = 3;
   int lda = 3;
   float A[] = { -58, -95.6, 44.3, 61.3, 50.3, 62.9, -41.9, 39.8, -65.6 };
   float X[] = { 73.6, -31.8, 33.1 };
   int incX = 1;
   float x_expected[] = { -3262.63, 1285.58, 33.1 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 230)");
     }
   };
  };

