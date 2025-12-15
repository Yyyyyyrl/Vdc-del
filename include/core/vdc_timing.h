#ifndef VDC_TIMING_H
#define VDC_TIMING_H

#include <string>
#include <chrono>
#include <map>
#include <vector>
#include <memory>

//! @brief Timer node for hierarchical timing measurement.
/*!
 * Each timer node represents a timing region in the program and can have child timers
 * to create a hierarchical timing tree.
 */
struct TimerNode {
    std::string name;                                           //!< Name of this timer
    std::chrono::high_resolution_clock::time_point start_time;  //!< Time point when timer started
    double elapsed_seconds;                                     //!< Total elapsed time in seconds
    bool is_running;                                            //!< Flag indicating if timer is currently running
    TimerNode* parent;                                          //!< Pointer to parent timer node
    std::vector<std::unique_ptr<TimerNode>> children;           //!< Child timer nodes

    //! @brief Constructor to create a new timer node.
    /*!
     * @param name Name of the timer
     * @param parent Pointer to parent timer node (nullptr for root)
     */
    TimerNode(const std::string& name, TimerNode* parent = nullptr);

    //! @brief Starts this timer.
    void start();

    //! @brief Stops this timer and accumulates elapsed time.
    void stop();

    //! @brief Gets the total elapsed time for this timer.
    /*!
     * @return Total elapsed time in seconds
     */
    double getElapsed() const;
};

//! @brief Singleton class for managing hierarchical timing measurements.
/*!
 * This class provides a singleton interface for timing various operations in the program.
 * Timers can be organized hierarchically to show parent-child relationships.
 */
class TimingStats {
public:
    //! @brief Gets the singleton instance of TimingStats.
    /*!
     * @return Reference to the singleton TimingStats instance
     */
    static TimingStats& getInstance();

    // Prevent copying and assignment
    TimingStats(const TimingStats&) = delete;
    TimingStats& operator=(const TimingStats&) = delete;

    //! @brief Starts a timer with the given name.
    /*!
     * If parent is specified, this timer becomes a child of that parent.
     * If parent is empty, it's a top-level timer.
     *
     * @param name Name of the timer to start
     * @param parent Name of the parent timer (empty string for top-level)
     */
    void startTimer(const std::string& name, const std::string& parent = "");

    //! @brief Stops the timer with the given name.
    /*!
     * @param name Name of the timer to stop
     * @param parent Name of the parent timer (empty string for top-level). Used to disambiguate identical child names under different parents.
     */
    void stopTimer(const std::string& name, const std::string& parent = "");

    //! @brief Prints the hierarchical timing report to stdout.
    void printReport() const;

    //! @brief Resets all timers.
    /*!
     * Clears all timer data. Useful for testing or resetting between runs.
     */
    void reset();

private:
    TimingStats();
    ~TimingStats() = default;

    struct TimerKey
    {
        TimerNode* parent;
        std::string name;

        bool operator<(const TimerKey& other) const
        {
            if (parent != other.parent) return parent < other.parent;
            return name < other.name;
        }
    };

    TimerNode* resolveParent(const std::string& parent);
    TimerNode* findTimer(const std::string& name);
    TimerNode* findTimerInSubtree(TimerNode* node, const std::string& name);
    void printNode(const TimerNode* node, int indent, bool is_last_child, const std::vector<bool>& ancestor_continues) const;
    std::string formatTime(double seconds) const;

    std::unique_ptr<TimerNode> root_;
    std::map<TimerKey, TimerNode*> timer_map_;  // For fast lookup by parent/name pair
};

//! @brief RAII helper for automatic timer start/stop.
/*!
 * This class automatically starts a timer on construction and stops it on destruction,
 * making it convenient to time scoped code blocks.
 */
class ScopedTimer {
public:
    //! @brief Constructor that starts a timer.
    /*!
     * @param name Name of the timer
     * @param parent Name of the parent timer (empty string for top-level)
     */
    ScopedTimer(const std::string& name, const std::string& parent = "");

    //! @brief Destructor that stops the timer.
    ~ScopedTimer();

private:
    std::string name_;  //!< Name of this scoped timer
    std::string parent_; //!< Name of the parent timer
};

#endif // VDC_TIMING_H
