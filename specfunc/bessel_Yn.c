float bessy(int n, float x)
Returns the Bessel function Y n (x) for positive x and n  2.
{
float bessy0(float x);
float bessy1(float x);
void nrerror(char error_text[]);
int j;
float by,bym,byp,tox;
if (n < 2) nrerror("Index n less than 2 in bessy");
tox=2.0/x;
by=bessy1(x); Starting values for the recurrence.
bym=bessy0(x);
for (j=1;j<n;j++) { Recurrence (6.5.7).
byp=j*tox*by-bym;
bym=by;
by=byp;
}
return by;
}
