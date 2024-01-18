 

    def get_allowed_number_of_cores_for_grid(self, grid, min_procs, max_procs):

        nx = grid.nx - 1
        ny = grid.ny - 1
        nz = grid.nz - 1

        anc = list()
        for np_used in range(min_procs, max_procs+1):

            for npex in range(1, np_used+1):

                allowed = True

                # strict matching in nodes (all nodes are used completely)
                if self.tasks_per_node != 0:
                    if np_used % self.tasks_per_node != 0:
                        allowed = False

                if np_used % npex != 0:
                    allowed = False

                npey = int(np_used / npex)

                # if fixed ratio is requested
                if self.npex_d_npey != 0.0 and npex / npey != self.npex_d_npey:
                    allowed = False

                # nx restlos teilbar auf npex
                if (nx + 1) % npex != 0:
                    allowed = False
                # ny restlos teilbar auf npey
                if (ny + 1) % npey != 0:
                    allowed = False

                if self.nx_d_ny != 0.0 and float(nx) / float(ny) != float(self.nx_d_ny):
                    allowed = False

                if self.multi_grid is True and (nx + 1) % 2 != 0:
                    allowed = False
                if self.multi_grid is True and (ny + 1) % 2 != 0:
                    allowed = False
                if self.multi_grid is True and nz % 2 != 0:
                    allowed = False

                if self.switch is True and (nx + 1) % npey != 0:
                    allowed = False
                if self.switch is True and (ny + 1) % npex != 0:
                    allowed = False

                if self.poisfft is True and nz % npex != 0:
                    allowed = False

                if npex > nx:
                    allowed = False
                if npey > ny:
                    allowed = False

                if self.spectra is True and nz % npey != 0:
                    allowed = False

                if self.tempterton is True and nx > 1 and ny > 1 and self.temperton_check(nx, ny) == 0:
                    allowed = False

                npxnpy = format(float(npex) / float(npey), '.2f')
                nxpex = float(nx + 1) / float(npex)
                nypey = float(ny + 1) / float(npey)

                if allowed:
                    anc.append((nx, ny, nz, npex, npey, npxnpy, (npex * npey), (nx * ny * nz), nxpex, nypey))
        return anc
