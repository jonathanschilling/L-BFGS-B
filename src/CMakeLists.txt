list(APPEND liblbfgsb_src
    ${CMAKE_CURRENT_SOURCE_DIR}/active.f
    ${CMAKE_CURRENT_SOURCE_DIR}/bmv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/cauchy.f
    ${CMAKE_CURRENT_SOURCE_DIR}/cmprlb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dcsrch.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dcstep.f
    ${CMAKE_CURRENT_SOURCE_DIR}/errclb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/formk.f
    ${CMAKE_CURRENT_SOURCE_DIR}/formt.f
    ${CMAKE_CURRENT_SOURCE_DIR}/freev.f
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsolb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/lnsrlb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/mainlb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/matupd.f
    ${CMAKE_CURRENT_SOURCE_DIR}/prn1lb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/prn2lb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/prn3lb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/projgr.f
    ${CMAKE_CURRENT_SOURCE_DIR}/setulb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/subsm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.f
)
set (liblbfgsb_src "${liblbfgsb_src}" PARENT_SCOPE)
