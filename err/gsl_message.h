#ifndef GSL_MESSAGE_H
#define GSL_MESSAGE_H

/* Provide a general messaging service for client use.  Messages can
 * be selectively turned off at compile time by defining an
 * appropriate message mask. Client code which uses the GSL_MESSAGE()
 * macro must provide a mask which is or'ed with the GSL_MESSAGE_MASK.
 *
 * The messaging service can be completely turned off
 * by defining GSL_MESSAGING_OFF.  */

void gsl_message(const char * message, const char * file, int line,
		 unsigned int mask);

#ifndef GSL_MESSAGE_MASK
#define GSL_MESSAGE_MASK 0xffffffffu /* default all messages allowed */
#endif

extern unsigned int gsl_message_mask ;

/* Provide some symolic masks for client ease of use. */

enum {
  GSL_MESSAGE_MASK_A = 1,
  GSL_MESSAGE_MASK_B = 2,
  GSL_MESSAGE_MASK_C = 4,
  GSL_MESSAGE_MASK_D = 8,
  GSL_MESSAGE_MASK_E = 16,
  GSL_MESSAGE_MASK_F = 32,
  GSL_MESSAGE_MASK_G = 64,
  GSL_MESSAGE_MASK_H = 128
} ;

#ifdef GSL_MESSAGING_OFF        /* throw away messages */ 
#define GSL_MESSAGE(message, mask) do { } while(0)
#else                           /* output all messages */
#define GSL_MESSAGE(message, mask) \
       do { \
       if (mask & GSL_MESSAGE_MASK) \
	 gsl_message (message, __FILE__, __LINE__, mask) ; \
       } while (0)
#endif

#endif /* GSL_MESSAGE_H */


