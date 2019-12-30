// @flow

import React from 'react';
import TableBody from '@material-ui/core/TableBody';
import TableRow from '@material-ui/core/TableRow';
import TableCell from '@material-ui/core/TableCell';
import type { ReadOnlyData, RowActionType, TableColumn } from './types';
import RowActions from './RowActions';

function Cell(
  column: TableColumn,
  row: ReadOnlyData,
  keyBase: string,
  actions: RowActionType[],
  keyField: string,
  size: string
) {
  if (column !== 'actions') {
    const value = row[column.id];
    return (
      <TableCell key={`${keyBase}-${column.id}`} align={column.align}>
        {column.format ? column.format(row) : value}
      </TableCell>
    );
  }
  return (
    <TableCell key={`${keyBase}-actions`} align="center">
      <RowActions actions={actions} data={row} keyField={keyField} size={size} />
    </TableCell>
  );
}

type Props = {
  isLoading: boolean,
  data: ReadOnlyData[],
  keyField: string,
  columns: TableColumn[],
  actions: RowActionType[],
  size: string
};

export default function Body({
  isLoading,
  data,
  keyField,
  columns,
  actions,
  size
}: Props) {
  return (
    <TableBody>
      {!isLoading &&
        data.map(row => (
          <TableRow
            hover
            role="checkbox"
            tabIndex={-1}
            key={`row-${row[keyField]}`}
          >
            {columns.map(column =>
              Cell(column, row, `row-${row[keyField]}`, actions, keyField, size)
            )}
          </TableRow>
        ))}
    </TableBody>
  );
}
