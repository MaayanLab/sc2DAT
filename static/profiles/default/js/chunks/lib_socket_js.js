"use strict";
(self["webpackChunkappyter_js"] = self["webpackChunkappyter_js"] || []).push([["lib_socket_js"],{

/***/ "./lib/socket.js":
/*!***********************!*\
  !*** ./lib/socket.js ***!
  \***********************/
/***/ (function(__unused_webpack___webpack_module__, __webpack_exports__, __webpack_require__) {

__webpack_require__.r(__webpack_exports__);
/* harmony import */ var uuid__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! uuid */ "./node_modules/uuid/dist/esm-browser/v4.js");
/* harmony import */ var socket_io_client__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! socket.io-client */ "./node_modules/socket.io-client/build/esm/index.js");


const socket = (0,socket_io_client__WEBPACK_IMPORTED_MODULE_0__["default"])({
  path: `${window._config.PREFIX}/socket.io/`,
  extraHeaders: {
    Authorization: (0,uuid__WEBPACK_IMPORTED_MODULE_1__["default"])()
  }
});
/* harmony default export */ __webpack_exports__["default"] = (socket);

/***/ })

}]);
//# sourceMappingURL=lib_socket_js.js.map