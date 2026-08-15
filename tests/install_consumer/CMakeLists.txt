cmake_minimum_required(VERSION 3.16)
project(rabbitsketch_install_consumer LANGUAGES CXX)

find_package(RabbitSketch 2 CONFIG REQUIRED)
add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE RabbitSketch::rabbitsketch_static)
add_executable(consumer_shared main.cpp)
target_link_libraries(consumer_shared PRIVATE RabbitSketch::rabbitsketch)
