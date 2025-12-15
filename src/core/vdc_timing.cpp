#include "core/vdc_timing.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

// ============================================================================
// TimerNode Implementation
// ============================================================================

TimerNode::TimerNode(const std::string& name, TimerNode* parent)
    : name(name), elapsed_seconds(0.0), is_running(false), parent(parent) {
}

void TimerNode::start() {
    if (!is_running) {
        start_time = std::chrono::high_resolution_clock::now();
        is_running = true;
    }
}

void TimerNode::stop() {
    if (is_running) {
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - start_time;
        elapsed_seconds += diff.count();
        is_running = false;
    }
}

double TimerNode::getElapsed() const {
    if (is_running) {
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = current_time - start_time;
        return elapsed_seconds + diff.count();
    }
    return elapsed_seconds;
}

// ============================================================================
// TimingStats Implementation
// ============================================================================

TimingStats::TimingStats() : root_(std::make_unique<TimerNode>("ROOT")) {
    timer_map_[{root_.get(), "ROOT"}] = root_.get();
}

TimingStats& TimingStats::getInstance() {
    static TimingStats instance;
    return instance;
}

void TimingStats::startTimer(const std::string& name, const std::string& parent) {
    TimerNode* parent_node = resolveParent(parent);
    TimerKey key{parent_node, name};

    // Check if timer already exists for this parent
    auto it = timer_map_.find(key);
    if (it != timer_map_.end()) {
        it->second->start();
        return;
    }

    // Create new timer
    auto new_timer = std::make_unique<TimerNode>(name, parent_node);
    TimerNode* timer_ptr = new_timer.get();
    parent_node->children.push_back(std::move(new_timer));
    timer_map_[key] = timer_ptr;
    timer_ptr->start();
}

void TimingStats::stopTimer(const std::string& name, const std::string& parent) {
    TimerNode* parent_node = resolveParent(parent);
    TimerKey key{parent_node, name};

    auto it = timer_map_.find(key);
    if (it != timer_map_.end()) {
        it->second->stop();
        return;
    }

    // Fallback: search by name if parent-specific timer was not found
    TimerNode* timer = findTimer(name);
    if (timer) timer->stop();
}

TimerNode* TimingStats::findTimer(const std::string& name) {
    return findTimerInSubtree(root_.get(), name);
}

TimerNode* TimingStats::findTimerInSubtree(TimerNode* node, const std::string& name) {
    if (node->name == name) {
        return node;
    }
    for (const auto& child : node->children) {
        TimerNode* result = findTimerInSubtree(child.get(), name);
        if (result) {
            return result;
        }
    }
    return nullptr;
}

std::string TimingStats::formatTime(double seconds) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << seconds;
    return oss.str();
}

void TimingStats::printNode(const TimerNode* node, int indent, bool is_last_child, const std::vector<bool>& ancestor_continues) const {
    if (node->name == "ROOT") {
        // Don't print the root node itself, just its children (the 8 main sections)
        std::vector<bool> empty_ancestors;
        for (size_t i = 0; i < node->children.size(); ++i) {
            printNode(node->children[i].get(), 0, i == node->children.size() - 1, empty_ancestors);
        }
        return;
    }

    // Calculate total time for children
    double children_total = 0.0;
    for (const auto& child : node->children) {
        children_total += child->getElapsed();
    }

    double elapsed = node->getElapsed();

    // Build the tree prefix
    std::string prefix = "[TIMING] ";

    if (indent > 0) {
        // Add one space for each indent level as base spacing
        for (int i = 0; i < indent; ++i) {
            prefix += " ";
        }

        // Add the tree character for this level
        if (is_last_child) {
            prefix += "└─ ";
        } else {
            prefix += "├─ ";
        }
    }

    // Print current node
    std::cout << prefix << node->name << " (" << formatTime(elapsed) << " s)";

    std::cout << "\n";

    // Print children with updated ancestor state
    if (!node->children.empty()) {
        std::vector<bool> new_ancestors = ancestor_continues;
        new_ancestors.push_back(!is_last_child);  // This level continues if we're not the last child

        for (size_t i = 0; i < node->children.size(); ++i) {
            printNode(node->children[i].get(), indent + 1, i == node->children.size() - 1, new_ancestors);
        }
    }
}

void TimingStats::printReport() const {
    std::cout << "\n====Execution Timing Stats====\n";

    std::vector<bool> empty_ancestors;

    // Find "Total Processing" node and print its children as top-level
    for (const auto& child : root_->children) {
        if (child->name == "Total Processing") {
            // Print the 8 main sections as top-level items
            for (size_t i = 0; i < child->children.size(); ++i) {
                printNode(child->children[i].get(), 0, i == child->children.size() - 1, empty_ancestors);
            }
            // Print total time from "Total Processing"
            std::cout << "==============================\n";
            std::cout << "Total Processing Time: " << formatTime(child->getElapsed()) << " seconds\n";
            std::cout << "==============================\n";
            return;
        }
    }

    // Fallback: print normally if "Total Processing" not found
    printNode(root_.get(), 0, true, empty_ancestors);
    double total = 0.0;
    for (const auto& child : root_->children) {
        total += child->getElapsed();
    }
    std::cout << "==============================\n";
    std::cout << "Total Processing Time: " << formatTime(total) << " seconds\n";
    std::cout << "==============================\n";
}

void TimingStats::reset() {
    timer_map_.clear();
    root_ = std::make_unique<TimerNode>("ROOT");
    timer_map_[{root_.get(), "ROOT"}] = root_.get();
}

TimerNode* TimingStats::resolveParent(const std::string& parent) {
    if (parent.empty()) {
        return root_.get();
    }
    TimerNode* parent_node = findTimer(parent);
    return parent_node ? parent_node : root_.get();
}

// ============================================================================
// ScopedTimer Implementation
// ============================================================================

ScopedTimer::ScopedTimer(const std::string& name, const std::string& parent)
    : name_(name), parent_(parent) {
    TimingStats::getInstance().startTimer(name_, parent_);
}

ScopedTimer::~ScopedTimer() {
    TimingStats::getInstance().stopTimer(name_, parent_);
}
