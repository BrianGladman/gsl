


void
  msg_checking (const char *test_description);

void
  msg_checking_params (const char *params_description,
		       const char *test_description);

void
  msg_result_status (int status);

void
  msg_result (const char *result_description);

void
  msg_error (const char *error_description);

void
  msg_summary (unsigned int tests, unsigned int passed, unsigned int failed);
