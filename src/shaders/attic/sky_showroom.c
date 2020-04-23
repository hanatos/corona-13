/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "spectrum.h"

extern float sky(rayhit_t *hit, float *dir, void *data)
{
  if(dir[1] < 0.0) return .0f;
  else
  {
    int d = (int)(dir[0]*10.0f/dir[1]);
    //if(abs(d) < 10 && d & 1)
    if(d < 5 && d > -5) return dir[1]*dir[1]*dir[1]*dir[1]*10.0f;
    else return 0.0f;
  }
}

