#!/usr/bin/env python
from PyQt5.QtCore import pyqtSignal, QThread
from math import floor
from palm_gf_files.palm_gf_conf import cfg_read
from datetime import datetime
import sqlite3


class MyThread(QThread):
    change_value = pyqtSignal(int)
    completed = pyqtSignal(int, int)
    busy = pyqtSignal()

    def run(self):
        # has to be defined in QThread class
        # starts calculation
        self.mainroutine()

    def mainroutine(self):
        # determines all possible solutions for given parameters
        start = datetime.now()
        conn = sqlite3.connect("palm_gf_data.db")
        c = conn.cursor()
        c.execute("DROP TABLE IF EXISTS grid_current")
        c.execute("DROP TABLE IF EXISTS grid_limits")
        c.execute("CREATE TABLE IF NOT EXISTS grid_current(nx INT, ny INT, nz INT, npex INT, npey INT, pxpy FLOAT,"
                  "np INT, ngpts FLOAT, nxnpex FLOAT, nynpey FLOAT)")
        c.execute("CREATE TABLE IF NOT EXISTS grid_limits(nx INT, ny INT, nz INT, npex INT, npey INT, pxpy FLOAT,"
                  "np INT, ngpts FLOAT, nxnpex FLOAT, nynpey FLOAT)")
        conn.commit()
        var = cfg_read()

        for j in range(12, 17):
            if var[j] == "True":
                var[j] = True
            else:
                var[j] = False

        nx_min = int(var[0])
        nx_max = int(var[1])
        ny_min = int(var[2])
        ny_max = int(var[3])
        nz_min = int(var[4])
        nz_max = int(var[5])
        procs_min = int(var[6])
        procs_max = int(var[7])
        tpn = int(var[8])
        dnxny = float(var[9])
        dpxpy = float(var[10])
        dpxpy_dev = int(var[11])
        poisfft = var[12]
        switch = var[13]
        temperton = var[14]
        mlt_grid = var[15]
        spectr = var[16]

        np_used = procs_min
        counter = 0
        cnt = 0
        nx = nx_min
        ny = ny_min
        nz = nz_min

        def factors(n):
            # prime factorization
            result = []
            for i in range(2, n + 1):
                s = 0
                while n / i == floor(n / float(i)):
                    n = n / float(i)
                    s += 1
                if s > 0:
                    for k in range(s):
                        result.append(i)
                        if n == 1:
                            return result

        def temperton_check(nx_, ny_):
            # uses previously defined prime factorization and determines if ny+1/nx+1 is multiple of 2,3,5
            nx_list = factors(nx_ + 1)
            ny_list = factors(ny_ + 1)
            return_value = 1
            i = 0
            if nx_list is not None:
                while i < len(nx_list):
                    if not (nx_list[i] == 2 or nx_list[i] == 3 or nx_list[i] == 5):
                        return_value = 0
                    i += 1
                j = 0
                if ny_list is not None:
                    while j < len(ny_list):
                        if not (ny_list[j] == 2 or ny_list[j] == 3 or ny_list[j] == 5):
                            return_value = 0
                        j += 1

            if nx_list is None:
                return_value = 0
            if ny_list is None:
                return_value = 0

            return return_value
            # * * *

        while np_used <= procs_max:
            a = 1
            while a <= np_used:
                if cnt != ((np_used - procs_min) / (procs_max - procs_min)) * 100:
                    cnt = ((np_used - procs_min) / (procs_max - procs_min)) * 100
                    self.change_value.emit(cnt)
                procs_var = np_used % a
                if procs_var != 0:
                    a += 1
                elif procs_var == 0:
                    npex = a
                    npey = int(np_used / npex)

                    if tpn != 0:
                        if np_used % tpn != 0:
                            a += 1
                            continue

                    if dpxpy != 0.:
                        if (dpxpy - dpxpy_dev*dpxpy) > (float(npex)/float(npey)) or \
                           (float(npex)/float(npey)) > (dpxpy + dpxpy_dev*dpxpy):
                            a += 1
                            continue

                    # if dpxpy != 0. and float(npex) / float(npey) != dpxpy:
                        # a += 1
                        # continue

                    while nx <= nx_max:
                        if (nx + 1) % npex != 0:
                            nx += 1
                            continue

                        if mlt_grid is True and nx % 2 != 0:
                            nx += 1
                            continue

                        if switch is True and (nx + 1) % npey != 0:
                            nx += 1
                            continue

                        if npex > nx:
                            nx += 1
                            continue

                        while ny <= ny_max:

                            if dnxny != 0. and float(nx) / float(ny) != float(dnxny):
                                ny += 1
                                continue

                            if (ny + 1) % npey != 0:
                                ny += 1
                                continue

                            if mlt_grid is True and ny % 2 != 0:
                                ny += 1
                                continue

                            if switch is True and (ny + 1) % npex != 0:
                                ny += 1
                                continue

                            if npey > ny:
                                ny += 1
                                continue

                            while nz <= nz_max:

                                if mlt_grid is True and nz % 2 != 0:
                                    nz += 1
                                    continue

                                if poisfft is True and nz % npex != 0:
                                    nz += 1
                                    continue

                                if spectr is True and nz % npey != 0:
                                    nz += 1
                                    continue

                                if temperton is True and nx > 1 and ny > 1 and temperton_check(nx, ny) == 0:
                                    nz += 1
                                    continue

                                npxnpy = format(float(npex) / float(npey), '.2f')
                                nxnpex = float(nx + 1) / float(npex)
                                nynpey = float(ny + 1) / float(npey)

                                c.execute("""INSERT OR REPLACE INTO grid_current(nx, ny, nz, npex, npey, pxpy,
                                np, ngpts, nxnpex, nynpey) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                                          (nx, ny, nz, npex, npey, npxnpy, (npex * npey), (nx * ny * nz),
                                           nxnpex, nynpey))

                                counter += 1
                                nz += 1
                            nz = nz_min
                            ny += 1
                        ny = ny_min
                        nx += 1
                    nx = nx_min
                    a += 1
            np_used += 1
            conn.commit()

        conn.commit()
        self.busy.emit()
        # read out limits of data in .db
        try:
            c.execute("SELECT nx FROM grid_current ORDER BY nx DESC LIMIT 1")
            mx_nx = c.fetchone()[0]
            c.execute("SELECT nx FROM grid_current ORDER BY nx  LIMIT 1")
            mn_nx = c.fetchone()[0]
            c.execute("SELECT ny FROM grid_current ORDER BY ny DESC LIMIT 1")
            mx_ny = c.fetchone()[0]
            c.execute("SELECT ny FROM grid_current ORDER BY ny  LIMIT 1")
            mn_ny = c.fetchone()[0]
            c.execute("SELECT nz FROM grid_current ORDER BY nz DESC LIMIT 1")
            mx_nz = c.fetchone()[0]
            c.execute("SELECT nz FROM grid_current ORDER BY nz  LIMIT 1")
            mn_nz = c.fetchone()[0]
            c.execute("SELECT npex FROM grid_current ORDER BY npex DESC LIMIT 1")
            mx_npex = c.fetchone()[0]
            c.execute("SELECT npex FROM grid_current ORDER BY npex  LIMIT 1")
            mn_npex = c.fetchone()[0]
            c.execute("SELECT npey FROM grid_current ORDER BY npey DESC LIMIT 1")
            mx_npey = c.fetchone()[0]
            c.execute("SELECT npey FROM grid_current ORDER BY npey  LIMIT 1")
            mn_npey = c.fetchone()[0]
            c.execute("SELECT pxpy FROM grid_current ORDER BY pxpy DESC LIMIT 1")
            mx_npxnpy = c.fetchone()[0]
            c.execute("SELECT pxpy FROM grid_current ORDER BY pxpy  LIMIT 1")
            mn_npxnpy = c.fetchone()[0]
            c.execute("SELECT np FROM grid_current ORDER BY np DESC LIMIT 1")
            mx_np = c.fetchone()[0]
            c.execute("SELECT np FROM grid_current ORDER BY np LIMIT 1")
            mn_np = c.fetchone()[0]
            c.execute("SELECT ngpts FROM grid_current ORDER BY ngpts DESC LIMIT 1")
            mx_ngpts = c.fetchone()[0]
            c.execute("SELECT ngpts FROM grid_current ORDER BY ngpts LIMIT 1")
            mn_ngpts = c.fetchone()[0]
            c.execute("SELECT nxnpex FROM grid_current ORDER BY nxnpex DESC LIMIT 1")
            mx_nxpex = c.fetchone()[0]
            c.execute("SELECT nxnpex FROM grid_current ORDER BY nxnpex LIMIT 1")
            mn_nxpex = c.fetchone()[0]
            c.execute("SELECT nynpey FROM grid_current ORDER BY nynpey DESC LIMIT 1")
            mx_nypey = c.fetchone()[0]
            c.execute("SELECT nynpey FROM grid_current ORDER BY nynpey LIMIT 1")
            mn_nypey = c.fetchone()[0]

            conn.commit()
            c.execute(
                """INSERT OR REPLACE INTO grid_limits(nx, ny, nz, npex, npey,
                 pxpy, np, ngpts, nxnpex, nynpey) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (mn_nx, mn_ny, mn_nz, mn_npex, mn_npey, mn_npxnpy, mn_np,  mn_ngpts, mn_nxpex, mn_nypey))

            c.execute(
                """INSERT OR REPLACE INTO grid_limits(nx, ny, nz, npex, npey,
                 pxpy, np, ngpts, nxnpex, nynpey) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (mx_nx, mx_ny, mx_nz, mx_npex, mx_npey, mx_npxnpy, mx_np, mx_ngpts, mx_nxpex, mx_nypey))
            conn.commit()

            c.close()
            conn.close()
        except TypeError:
            pass

        finish = datetime.now()
        timedelta = finish - start
        timer = (int(timedelta.seconds))
        self.completed.emit(counter, timer)
