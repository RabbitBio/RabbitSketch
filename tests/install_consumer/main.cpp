#include "api/RuntimeInfo.h"
#include "api/Version.h"
#include "estimators/UnifiedQuery.h"
#include "fastkmv.h"

#include <iostream>

int main() {
    Sketch::FastKMV left(64, 5, 42);
    Sketch::FastKMV right(32, 5, 42);
    const char sequence[] = "ACGTACGTACGT";
    left.update(sequence, sizeof(sequence) - 1);
    right.update(sequence, sizeof(sequence) - 1);
    const auto result = Sketch::Query::query(left, right);
    std::cout << Sketch::API::VERSION << ' '
              << Sketch::Runtime::runtimeInfo().selected_byte_path << ' '
              << left.size() << ' ' << result.value << '\n';
    return left.size() == 0 || !result.estimate_available ||
        result.value != 1.0 ? 1 : 0;
}
