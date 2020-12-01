/*
 * Copyright (C) 2020 Stefan Westerfeld
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef AUDIOWMARK_THREAD_POOL_HH
#define AUDIOWMARK_THREAD_POOL_HH

#include <vector>
#include <thread>
#include <functional>
#include <set>
#include <mutex>
#include <condition_variable>

class ThreadPool
{
  std::vector<std::thread> threads;

  struct Job
  {
    std::function<void()> fun;
  };

  std::mutex                mutex;
  std::condition_variable   cond;
  std::condition_variable   main_cond;
  std::vector<Job>          jobs;
  size_t                    jobs_added = 0;
  size_t                    jobs_done = 0;
  bool                      stop_workers = false;

  bool worker_next_job (Job& job);
  void worker_run();

public:
  ThreadPool();
  ~ThreadPool();

  void add_job (std::function<void()> fun);
  void wait_all();
};

#endif /* AUDIOWMARK_THREAD_POOL_HH */
