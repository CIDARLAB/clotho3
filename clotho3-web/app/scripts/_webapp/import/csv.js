angular.module('clotho.webapp')
.service('ImportCSVService', function ($window, $q, ClothoExtensions) {

  var self = this;

  // baby parse (https://github.com/Rich-Harris/BabyParse) is a fork of PapaParse (https://github.com/mholt/PapaParse)
  // it is smaller, but cannot use:
  //    downloading (do yourself)
  //    service workers (may be an issue for large files)
  //    encoding
  //    chunk

  self.ready = ClothoExtensions.mixin('lib/babyparse.js')
  .then(function () {
    console.log($window.Baby);
    return $window.Baby;
  });

  self.defaults = {header: true};

    /**
     * @name parse
     * @description
     * Parse a CSV file to a JSON object using BabyParse
     *
     * @param csv {String} CSV to parse
     * @param options {Object} params for parsing
     * @returns {Object} object with fields data, errors, meta
     */
  self.parse = function (csv, options) {

    if (!angular.isString(csv)) {
      return $q.when({});
    }

    options = angular.extend({}, self.defaults, options);

    return self.ready.then(function (parser) {
      return parser.parse(csv, options);
    });
  };
})
.controller('ImportCSVCtrl', function ($scope, $q, Clotho, ImportCSVService) {

  $scope.csvText = '';
  $scope.options = {
    header: true
  };

  function parseCsv () {
    ImportCSVService.parse.apply(null, arguments)
    .then(function (parsed) {
      //note - actual data is in parsed.data
      $scope.parsed = parsed;
      $scope.parsedData = parsed.data;

      //set up field map
      //take first entry as representative
      if (parsed.data.length) {
        $scope.fieldDemo = parsed.data[0];
        $scope.fieldMap = {};
        angular.forEach(Object.keys($scope.fieldDemo), function (key) {
          $scope.fieldMap[key] = key;
        });
      }
    });
  }

  $scope.process = function () {
    if (angular.isDefined($scope.csvText) && $scope.csvText.length) {
      parseCsv($scope.csvText, $scope.options);
    }
  };

  $scope.removeField = function (fieldName) {
    delete $scope.fieldMap[fieldName];
  };

  //given array of objects and dictionary of old keys to new keys, create new object with new keys
  $scope.reparse = function () {
    var mapped = [];
    angular.forEach($scope.parsedData, function (item) {
      var temp = {};
      //todo - ensure only goes over own properties
      angular.forEach($scope.fieldMap, function (newkey, oldkey) {
        temp[newkey] = item[oldkey];
      });
      mapped.push(temp);
    });
    $scope.parsedData = mapped;
  };

  $scope.importParsed = function () {
    $scope.createdIds = [];
    angular.forEach($scope.parsedData, function (item) {
      Clotho.create(item).then(function (createdId) {
        $scope.createdIds.push(createdId);
      });
    });
  }
});
