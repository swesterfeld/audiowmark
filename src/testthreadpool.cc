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

#include <stdio.h>
#include <unistd.h>

#include "threadpool.hh"

int
main()
{
  ThreadPool tp;

  int result1 = 0;
  int result2 = 0;
  tp.add_job ([&result1](){printf ("A\n"); sleep (2); printf ("A done\n"); result1 = 123;});
  tp.add_job ([&result2](){printf ("B\n"); sleep (3); printf ("B done\n"); result2 = 456;});

  tp.wait_all();
  printf ("===\n");
  printf ("results: %d, %d\n", result1, result2);
}
