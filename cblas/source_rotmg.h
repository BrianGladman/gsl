/* blas/source_rotmg.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

{
  const BASE G = 4096.0, G2 = G * G;
  BASE D1 = *d1, D2 = *d2, x1 = *b1, y1 = b2;
  BASE h11, h12, h21, h22, u;

  BASE c, s;

  /* case of d1 < 0, appendix A, second to last paragraph */

  if (D1 < 0.0) {
    P[0] = -1;
    P[1] = 0;
    P[2] = 0;
    P[3] = 0;
    P[4] = 0;
    *d1 = 0;
    *d2 = 0;
    *b1 = 0;
    return;
  }

  if (D2 * y1 == 0.0) {
    P[0] = -2;			/* case of H = I, page 315 */
    return;
  }

  c = fabs(D1 * x1 * x1);
  s = fabs(D2 * y1 * y1);

  if (c > s) {
    /* case of equation A6 */

    P[0] = 0.0;

    h11 = 1;
    h12 = (D2 * y1) / (D1 * x1);
    h21 = -y1 / x1;
    h22 = 1;

    u = 1 - h21 * h12;

    if (u <= 0.0) {		/* the case u <= 0 is rejected */
      P[0] = -1;
      P[1] = 0;
      P[2] = 0;
      P[3] = 0;
      P[4] = 0;
      *d1 = 0;
      *d2 = 0;
      *b1 = 0;
      return;
    }

    D1 /= u;
    D2 /= u;
    x1 *= u;
  } else {
    /* case of equation A7 */

    if (D2 * y1 * y1 < 0.0) {
      P[0] = -1;
      P[1] = 0;
      P[2] = 0;
      P[3] = 0;
      P[4] = 0;
      *d1 = 0;
      *d2 = 0;
      *b1 = 0;
      return;
    }

    P[0] = 1;

    h11 = (D1 * x1) / (D2 * y1);
    h12 = 1;
    h21 = -1;
    h22 = x1 / y1;

    u = 1 + h11 * h22;

    D1 /= u;
    D2 /= u;

    {
      BASE tmp = D2;
      D2 = D1;
      D1 = tmp;
    }

    x1 = y1 * u;
  }

  /* rescale D1 to range [1/G2,G2] */

  while (D1 <= 1.0 / G2 && D1 != 0.0) {
    P[0] = -1;
    D1 *= G2;
    x1 /= G;
    h11 /= G;
    h12 /= G;
  }

  while (D1 >= G2) {
    P[0] = -1;
    D1 /= G2;
    x1 *= G;
    h11 *= G;
    h12 *= G;
  }

  /* rescale D2 to range [1/G2,G2] */

  while (fabs(D2) <= 1.0 / G2 && D2 != 0.0) {
    P[0] = -1;
    D2 *= G2;
    h21 /= G;
    h22 /= G;
  }

  while (fabs(D2) >= G2) {
    P[0] = -1;
    D2 /= G2;
    h21 *= G;
    h22 *= G;
  }

  *d1 = D1;
  *d2 = D2;
  *b1 = x1;

  if (P[0] == -1.0) {
    P[1] = h11;
    P[2] = h21;
    P[3] = h12;
    P[4] = h22;
  } else if (P[0] == 0.0) {
    P[2] = h21;
    P[3] = h12;
  } else if (P[0] == 1.0) {
    P[1] = h11;
    P[4] = h22;
  }
}
