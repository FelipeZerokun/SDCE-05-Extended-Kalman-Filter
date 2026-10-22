# Synthetic inputs keep CLI verification independent of the ignored Udacity data.
file(MAKE_DIRECTORY "${TEST_DIR}")
set(valid "${TEST_DIR}/stationary.txt")
file(WRITE "${valid}" "L 2 0 0 2 0 0 0 0 0\nR 2 0 0 100000 2 0 0 0 0 0\nL 2 0 200000 2 0 0 0 0 0\n")
execute_process(COMMAND "${REPLAY}" "${valid}" --check-rmse
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
string(FIND "${output}" "Measurements: 3 (lidar: 2, radar: 1)" counts_position)
if(NOT result EQUAL 0 OR counts_position EQUAL -1
        OR NOT output MATCHES "RMSE \\[px py vx vy\\]: +0\\.000000 +0\\.000000 +0\\.000000 +0\\.000000"
        OR NOT output MATCHES "PASS")
    message(FATAL_ERROR "Valid replay failed (${result}): ${output}${error}")
endif()

set(invalid "${TEST_DIR}/malformed.txt")
file(WRITE "${invalid}" "L 2 0 0 2 0 0 0 0 0\ninvalid\n")
execute_process(COMMAND "${REPLAY}" "${invalid}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 1 OR NOT error MATCHES ":2:")
    message(FATAL_ERROR "Malformed row must fail with its line number: ${result}, ${error}")
endif()

set(empty "${TEST_DIR}/empty.txt")
file(WRITE "${empty}" "")
execute_process(COMMAND "${REPLAY}" "${empty}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 1 OR NOT error MATCHES "no measurements")
    message(FATAL_ERROR "Empty input must fail: ${result}, ${error}")
endif()

set(inaccurate "${TEST_DIR}/inaccurate.txt")
file(WRITE "${inaccurate}" "L 2 0 0 20 0 0 0 0 0\n")
execute_process(COMMAND "${REPLAY}" "${inaccurate}" --check-rmse
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 2 OR NOT output MATCHES "FAIL")
    message(FATAL_ERROR "Failed benchmark must return 2: ${result}, ${output}${error}")
endif()
