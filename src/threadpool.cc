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

#include "threadpool.hh"
#include "utils.hh"

bool
ThreadPool::worker_next_job (Job& job)
{
  std::unique_lock<std::mutex> lck (mutex);

  if (jobs.empty())
    cond.wait (lck);

  if (jobs.empty())
    return false;

  job = jobs.front();
  jobs.erase (jobs.begin());
  return true;
}

void
ThreadPool::worker_run()
{
  while (!stop_workers)
    {
      Job job;
      if (worker_next_job (job))
        {
          job.fun();

          std::lock_guard<std::mutex> lg (mutex);
          jobs_done.insert (job.id);

          main_cond.notify_one();
        }
    }
}

ThreadPool::ThreadPool()
{
  for (unsigned int i = 0; i < std::thread::hardware_concurrency(); i++)
    {
      threads.push_back (std::thread (&ThreadPool::worker_run, this));
    }
}

int
ThreadPool::add_job (std::function<void()> fun)
{
  std::lock_guard<std::mutex> lg (mutex);
  Job job;
  job.fun = fun;
  job.id  = next_job_id++;
  jobs.push_back (job);

  cond.notify_one();

  return job.id;
}

void
ThreadPool::wait_jobs (std::vector<int>& ids)
{
  for (;;)
    {
      std::unique_lock<std::mutex> lck (mutex);

      /* check if at least one of the jobs is still running */
      bool done = true;
      for (auto id : ids)
        if (jobs_done.count (id) == 0)
          done = false;

      if (done)
        {
          for (auto id : ids)
            jobs_done.erase (id);
          return;
        }

      main_cond.wait (lck);
    }
}

ThreadPool::~ThreadPool()
{
  {
    std::lock_guard<std::mutex> lg (mutex);
    stop_workers = true;
    cond.notify_all();
  }

  for (auto& t : threads)
    t.join();

  if (jobs_done.size())
    {
      // user must wait for each job before deleting the ThreadPool
      error ("audiowmark: %zd open jobs in ThreadPool::~ThreadPool() - this should not happen\n", jobs_done.size());
    }
}
