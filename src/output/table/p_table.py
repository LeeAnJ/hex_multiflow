"""Style options
PrettyTable has a number of style options which control various aspects of how tables are displayed. You have the
freedom to set each of these options individually to whatever you prefer. The set_style method just does this
automatically for you.

The options are these:

border - A boolean option (must be True or False). Controls whether a border is drawn inside and around the table.
preserve_internal_border - A boolean option (must be True or False). Controls whether borders are still drawn within the table even when border=False.
header - A boolean option (must be True or False). Controls whether the first row of the table is a header showing the names of all the fields.
hrules - Controls printing of horizontal rules after rows. Allowed values: FRAME, HEADER, ALL, NONE - note that these are variables defined inside the prettytable module so make sure you import them or use prettytable.FRAME etc.
vrules - Controls printing of vertical rules between columns. Allowed values: FRAME, ALL, NONE.
int_format - A string which controls the way integer data is printed. This works like: print("%<int_format>d" % data)
float_format - A string which controls the way floating point data is printed. This works like: print("%<float_format>f" % data)
custom_format - A Dictionary of field and callable. This allows you to set any format you want pf.custom_format["my_col_int"] = ()lambda f, v: f"{v:,}". The type of the callable if callable[[str, Any], str]
padding_width - Number of spaces on either side of column data (only used if left and right paddings are None).
left_padding_width - Number of spaces on left-hand side of column data.
right_padding_width - Number of spaces on right-hand side of column data.
vertical_char - Single character string used to draw vertical lines. Default is |.
horizontal_char - Single character string used to draw horizontal lines. Default is -.
_horizontal_align_char - single character string used to indicate column alignment in horizontal lines. Default is : for Markdown, otherwise None.
junction_char - Single character string used to draw line junctions. Default is +.
top_junction_char - single character string used to draw top line junctions. Default is junction_char.
bottom_junction_char - single character string used to draw bottom line junctions. Default is junction_char.
right_junction_char - single character string used to draw right line junctions. Default is junction_char.
left_junction_char - single character string used to draw left line junctions. Default is junction_char.
top_right_junction_char - single character string used to draw top-right line junctions. Default is junction_char.
top_left_junction_char - single character string used to draw top-left line junctions. Default is junction_char.
bottom_right_junction_char - single character string used to draw bottom-right line junctions. Default is junction_char
bottom_left_junction_char - single character string used to draw bottom-left line junctions. Default is junction_char."""
import sys

from prettytable import PrettyTable
from dataclasses import dataclass
from prettytable import MSWORD_FRIENDLY, DEFAULT, PLAIN_COLUMNS, MARKDOWN, ORGMODE, SINGLE_BORDER, DOUBLE_BORDER

from src.hex.hex import Hex


def print_input_data_as_table(hex: Hex):
    table_input_data_header = Table.header_input_data(n_hot_flows=hex.pack.hot.n_flows,
                                                      n_cold_flows=hex.pack.cold.n_flows)
    table_input_data_data = Table.data_input_data(hex=hex,
                                                  n_hot_flows=hex.pack.hot.n_flows,
                                                  n_cold_flows=hex.pack.cold.n_flows)
    table_input_data = Table(data=TableInputData(title=' Table: Input data ',
                                                 header=table_input_data_header,
                                                 preheader=None,
                                                 columns=table_input_data_data))
    table_input_data.print()

def print_dt_with_cross_section_as_table(hex: Hex):
    table_title = f' Table: dt = f(i_cross_section) of ({hex.pack.hot.n_flows} hot x {hex.pack.cold.n_flows} cold) flows HEX with Eps = {hex.eps_hex}'
    table_preheader = None
    table_header = Table.header_dt_with_cross_section(n_hot_flows=hex.pack.hot.n_flows, n_cold_flows=hex.pack.cold.n_flows)
    table_data = Table.data_dt_with_cross_section(hex=hex)
    table_dt = Table(data=TableInputData(title=table_title, header=table_header, preheader=table_preheader, columns=table_data))
    table_dt.print()

    print(f'  dt_min inter HEXs packs,    [K] = {hex.dt_min_inter_flows[-1][0]: 7.4f}')
    print(f'  dt_max inter HEXs packs,    [K] = {hex.dt_max_inter_flows[-1][0]: 7.4f}')
    print(f'  dt_avr inter HEXs packs,    [K] = {hex.dt_avr_inter_flows[-1][0]: 7.4f}')

def print_t_dq_with_cross_section_as_table(hex: Hex, qt_is_operable_flag):
    table_title = f' Table: t, dq = f(i_cross_section) of ({hex.pack.hot.n_flows} hot x {hex.pack.cold.n_flows} cold) flows HEX with Eps = {hex.eps_hex}'
    table_preheader = Table.preheader_t_dq_with_cross_section(hex=hex)
    table_header = Table.header_t_dq_with_cross_section(n_hot_flows=hex.pack.hot.n_flows,
                                                        n_cold_flows=hex.pack.cold.n_flows)
    table_data = Table.data_t_dq_with_cross_section(hex=hex)
    table_qt = Table(
        data=TableInputData(title=table_title, preheader=table_preheader, header=table_header, columns=table_data))
    table_qt.print()
    print(f'  qt_diagram is operable:            {qt_is_operable_flag}')
    print(f'  Heat transferred  in HEX,   [W] = {hex.q_act_w: 7.4f}')
    print(f'  Entropy generated in HEX, [W/K] = {hex.ds_wk_sum: 7.4f}')

def print_t_ds_with_cross_section_as_table(hex: Hex):
    table_title = f' Table: t, ds = f(i_cross_section) of ({hex.pack.hot.n_flows} hot x {hex.pack.cold.n_flows} cold) flows HEX with Eps = {hex.eps_hex}'
    table_preheader = Table.preheader_t_ds_with_cross_section(hex=hex)
    table_header = Table.header_t_ds_with_cross_section(n_hot_flows=hex.pack.hot.n_flows,
                                                        n_cold_flows=hex.pack.cold.n_flows)
    table_data = Table.data_t_ds_with_cross_section(hex=hex)
    table_ts = Table(
        data=TableInputData(title=table_title, preheader=table_preheader, header=table_header, columns=table_data))
    table_ts.print()


@dataclass
class TableInputData:
    title: str
    preheader: str | None
    header: list
    columns: list


class Table:
    def __init__(self, data: TableInputData):
        self.table = PrettyTable(padding_width=1)
        self.table_title = data.title
        self.table.preheader = data.preheader
        self.table.field_names = data.header
        n_columns = len(self.table.field_names)

        self.table_data = list()
        if n_columns != len(data.columns):
            sys.exit('n_data_rows != n_header_fields in Table')

        n_rows = len(data.columns[0])
        for i in range(n_rows):
            row = list()
            if not isinstance(data.columns[0][i], str):
                row.append('{:.0f}'.format(data.columns[0][i]))
            else:
                row.append(data.columns[0][i])
            for j in range(1, n_columns):
                if not isinstance(data.columns[j][i], str):
                    row.append('{:.4f}'.format(data.columns[j][i]))
                else:
                    row.append(data.columns[j][i])
            self.table_data.append(row)

        self.table.add_rows(self.table_data)
        self.table.align = "c"

    def print(self):
        # self.table.align = "c"
        # self.table.align["dq_cold_sum, [W]"] = "r"

        # self.table.set_style(PLAIN_COLUMNS)         # хорошо для копирования данных в OriginLab
        # self.table.set_style(MSWORD_FRIENDLY)
        # self.table.set_style(DEFAULT)
        # self.table.set_style(MARKDOWN)
        # self.table.set_style(ORGMODE)
        self.table.set_style(SINGLE_BORDER)  # красиво выглядит, но копир. в OriginLab не выходит,
        # self.table.set_style(DOUBLE_BORDER)        # а в Notepad++ копируется хорошо

        print('')
        print(self.table_title)
        if self.table.preheader is None:
            print(self.table)
        else:
            print(self.table.get_string(title=self.table.preheader))

    # ------------------------------------------------------------------------------------------------------------------
    # Input data table
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def header_input_data(n_hot_flows, n_cold_flows):
        header = ["pack", "flow[i]", "component name", "composition, [kg/kg]/[mol/mol]",
                  "t_in, [K]/[oC]", "p, [kPa]/[bar]", "m, [kg/s]/[mol/s]"]
        return header

    @staticmethod
    def data_input_data(hex, n_hot_flows, n_cold_flows):
        colon_0 = []
        for i in range(n_hot_flows):
            colon_0.append('hot')
        for i in range(n_cold_flows):
            colon_0.append('cold')
        data = [colon_0]

        colon_1 = []
        for i in range(n_hot_flows):
            colon_1.append(f'flow[{i}]')
        for i in range(n_cold_flows):
            colon_1.append(f'flow[{i}]')
        data.append(colon_1)

        colon_2 = []
        for i in range(n_hot_flows):
            colon_2.append(hex.pack.hot.flow[i].fluid.get_name_compact())
        for i in range(n_cold_flows):
            colon_2.append(hex.pack.cold.flow[i].fluid.get_name_compact())
        data.append(colon_2)

        colon_3 = []
        for i in range(n_hot_flows):
            colon_3.append(f'({hex.pack.hot.flow[i].fluid.get_composition_kgkg_compact_rounded()}) / '
                           f'({hex.pack.hot.flow[i].fluid.get_composition_molmol_compact_rounded()})')
        for i in range(n_cold_flows):
            colon_3.append(f'({hex.pack.cold.flow[i].fluid.get_composition_kgkg_compact_rounded()}) / '
                           f'({hex.pack.cold.flow[i].fluid.get_composition_molmol_compact_rounded()})')

        data.append(colon_3)

        colon_4 = []
        for i in range(n_hot_flows):
            colon_4.append(f'{hex.pack.hot.flow[i].inlet.t_k} / {round(hex.pack.hot.flow[i].inlet.t_k - 273.15, 4)}')
        for i in range(n_cold_flows):
            colon_4.append(f'{hex.pack.cold.flow[i].inlet.t_k} / {round(hex.pack.cold.flow[i].inlet.t_k - 273.15, 4)}')
        data.append(colon_4)

        colon_5 = []
        for i in range(n_hot_flows):
            colon_5.append(f'{hex.pack.hot.flow[i].inlet.p_kpa} / {round(hex.pack.hot.flow[i].inlet.p_kpa / 100.0, 4)}')
        for i in range(n_cold_flows):
            colon_5.append(
                f'{hex.pack.cold.flow[i].inlet.p_kpa} / {round(hex.pack.cold.flow[i].inlet.p_kpa / 100.0, 4)}')
        data.append(colon_5)

        colon_6 = []
        for i in range(n_hot_flows):
            colon_6.append(f'{hex.pack.hot.flow[i].m_kgs} / {round(hex.pack.hot.flow[i].m_mols, 4)}')
        for i in range(n_cold_flows):
            colon_6.append(f'{hex.pack.cold.flow[i].m_kgs} / {round(hex.pack.cold.flow[i].m_mols, 4)}')
        data.append(colon_6)

        return data

    # ------------------------------------------------------------------------------------------------------------------
    # t, dq = f(i_cross_section) table
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def preheader_t_dq_with_cross_section(hex):
        if hex.pack.hot.has_min_q_w:  # hot pack is a "passive" one
            # cold packet: active
            eps = tuple([round(abs(hex.pack.cold.flow[i].dq_w[-1] / hex.pack.cold.dq_w[-1]), 2)
                         for i in range(hex.pack.cold.n_flows)])
            # hot packet: passive
            preheader = f'HOT packet of {hex.pack.hot.n_flows}-flows is "passive" and has INVARIABLE eps_set = ' \
                        f'{hex.pack.hot.eps_sets[0]} > ' \
                        f'COLD packet of {hex.pack.cold.n_flows}-flows is "active" and has CURRENT eps_set = {eps}'
        else:
            # hot packet: active
            eps = tuple([round(abs(hex.pack.hot.flow[i].dq_w[0] / hex.pack.hot.dq_w[0]), 2)
                         for i in range(hex.pack.hot.n_flows)])
            preheader = f'HOT of {hex.pack.hot.n_flows}-flows packet is "active" and has CURRENT eps_set = {eps} > ' \
                        f'COLD packet of {hex.pack.cold.n_flows}-flows is "passive" and has ' \
                        f'INVARIABLE eps_set = {hex.pack.cold.eps_sets[0]}'
        return preheader

    @staticmethod
    def header_t_dq_with_cross_section(n_hot_flows, n_cold_flows):
        # table_header = ["N", "t,[K]", "dq,[W]", "dq_sum,[W]", "t,[K]-", "dq,[W]-", "t,[K]--", "dq,[W]--", "dq_sum,[W]-"]
        header = ["CS[i]"]
        for i in range(n_hot_flows):
            header.append(f'T_hot[{i}], [K]')
            header.append(f'dQ_hot[{i}], [W]')
        header.append(f'dQ_hot_sum, [W]')
        header.append('>')
        for i in range(n_cold_flows):
            header.append(f'T_cold[{i}], [K]')
            header.append(f'dQ_cold[{i}], [W]')
        header.append(f'dQ_cold_sum, [W]')
        return header

    @staticmethod
    def data_t_dq_with_cross_section(hex):
        empty_cells = ['>'] * (hex.n_qt_nods)
        data = [[i for i in range(hex.n_qt_nods)]]
        for i in range(hex.pack.hot.n_flows):
            data.append(hex.pack.hot.flow[i].t_k)
            data.append(hex.pack.hot.flow[i].dq_w)
        data.append(hex.pack.hot.dq_w)
        data.append(empty_cells),
        for i in range(hex.pack.cold.n_flows):
            data.append(hex.pack.cold.flow[i].t_k)
            data.append(hex.pack.cold.flow[i].dq_w)
        data.append(hex.pack.cold.dq_w)
        return data

    # ------------------------------------------------------------------------------------------------------------------
    # dt = f(i_cross_section) table
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def header_dt_with_cross_section(n_hot_flows, n_cold_flows):
        # table_header = ["N", "t,[K]", "dq,[W]", "dq_sum,[W]", "t,[K]-", "dq,[W]-", "t,[K]--", "dq,[W]--", "dq_sum,[W]-"]
        header = ["CS[i]"]
        for i in range(n_hot_flows):
            for j in range(n_cold_flows):
                header.append(f'T_hot[{i}]-T_cold[{j}], [K]')
        return header

    @staticmethod
    def data_dt_with_cross_section(hex):
        # данные в таблицу добавляем сформированными столбцами
        empty_cells = ['>'] * (hex.n_qt_nods)
        index = [str(i) for i in range(hex.n_qt_nods)]
        index.append('------')
        index.append('min ->')
        index.append('max ->')
        index.append('avr ->')
        data = [index]  # индекс поперечного сечения потока (=теплообменника)
        for i in range(hex.pack.hot.n_flows):
            for j in range(hex.pack.cold.n_flows):
                column = list(hex.dt_inter_flows[i][j])
                column.append('---------------------')
                column.append(hex.dt_min_inter_flows[i][j])
                column.append(hex.dt_max_inter_flows[i][j])
                column.append(hex.dt_avr_inter_flows[i][j])
                data.append(column)

        return data

    # ------------------------------------------------------------------------------------------------------------------
    # t, ds = f(i_cross_section) table
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def preheader_t_ds_with_cross_section(hex):
        if hex.pack.hot.has_min_q_w:  # hot pack is a "passive" one
            # cold packet: active
            eps = tuple([round(abs(hex.pack.cold.flow[i].dq_w[-1] / hex.pack.cold.dq_w[-1]), 2)
                         for i in range(hex.pack.cold.n_flows)])
            # hot packet: passive
            preheader = f'HOT packet of {hex.pack.hot.n_flows}-flows is "passive" and has INVARIABLE eps_set = ' \
                        f'{hex.pack.hot.eps_sets[0]} > ' \
                        f'COLD packet of {hex.pack.cold.n_flows}-flows is "active" and has CURRENT eps_set = {eps}'
        else:
            # hot packet: active
            eps = tuple([round(abs(hex.pack.hot.flow[i].dq_w[0] / hex.pack.hot.dq_w[0]), 2)
                         for i in range(hex.pack.hot.n_flows)])
            preheader = f'HOT of {hex.pack.hot.n_flows}-flows packet is "active" and has CURRENT eps_set = {eps} > ' \
                        f'COLD packet of {hex.pack.cold.n_flows}-flows is "passive" and has ' \
                        f'INVARIABLE eps_set = {hex.pack.cold.eps_sets[0]}'
        return preheader

    @staticmethod
    def header_t_ds_with_cross_section(n_hot_flows, n_cold_flows):
        # table_header = ["N", "t,[K]", "dq,[W]", "dq_sum,[W]", "t,[K]-", "dq,[W]-", "t,[K]--", "dq,[W]--", "dq_sum,[W]-"]
        header = ["CS[i]"]
        for i in range(n_hot_flows):
            header.append(f'T_hot[{i}], [K]')
            header.append(f'dS_hot[{i}], [W/K]')
        header.append(f'dS_hot_sum, [W/K]')
        header.append('>')
        for i in range(n_cold_flows):
            header.append(f'T_cold[{i}], [K]')
            header.append(f'dS_cold[{i}], [W/K]')
        header.append(f'dS_cold_sum, [W/K]')
        header.append(f'dS_hex, [W/K]')
        return header

    @staticmethod
    def data_t_ds_with_cross_section(hex):
        empty_cells = ['>'] * (hex.n_qt_nods)
        data = [[i for i in range(hex.n_qt_nods)]]
        for i in range(hex.pack.hot.n_flows):
            data.append(hex.pack.hot.flow[i].t_k)
            data.append(hex.pack.hot.flow[i].ds_wk)
        data.append(hex.pack.hot.ds_wk)
        data.append(empty_cells),
        for i in range(hex.pack.cold.n_flows):
            data.append(hex.pack.cold.flow[i].t_k)
            data.append(hex.pack.cold.flow[i].ds_wk)
        data.append(hex.pack.cold.ds_wk)
        data.append(hex.ds_wk)
        return data
