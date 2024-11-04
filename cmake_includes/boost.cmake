if(BUILD_SHARED_LIBS)
    add_definitions(-DBOOST_ALL_DYN_LINK)
    add_definitions(-DBOOST_LOG_DYN_LINK)
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DBOOST_STACKTRACE_DYN_LINK)
    add_definitions(-DBOOST_STACKTRACE_LINK)
else()
    set(Boost_USE_STATIC_LIBS ON)
    # set(Boost_USE_STATIC_RUNTIME ON)
endif()
find_package(
    Boost REQUIRED
    COMPONENTS filesystem regex program_options thread log_setup log
    # signals2 is header-only lockfree is header-only
    OPTIONAL_COMPONENTS unit_test_framework timer stacktrace_basic stacktrace_backtrace stacktrace_windbg)
include_directories(${Boost_INCLUDE_DIRS})

set(ART_MODERN_LINK_LIBS
    ${ART_MODERN_LINK_LIBS}
    Boost::filesystem
    Boost::regex
    Boost::program_options
    Boost::thread
    Boost::log_setup
    Boost::log)
unset(WITH_BOOST_TIMER)
if(Boost_timer_FOUND)
    set(WITH_BOOST_TIMER ON)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::timer)
endif()
if(Boost_stacktrace_backtrace_FOUND)
    add_definitions(-DBOOST_STACKTRACE_USE_BACKTRACE)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::stacktrace_backtrace)
elseif(Boost_stacktrace_windbg_FOUND)
    add_definitions(-DBOOST_STACKTRACE_USE_WINDBG)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::stacktrace_windbg)
elseif(Boost_stacktrace_basic_FOUND)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::stacktrace_basic)
endif()
