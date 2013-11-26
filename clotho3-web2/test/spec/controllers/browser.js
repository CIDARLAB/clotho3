'use strict';

describe('Controller: BrowserCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var BrowserCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    BrowserCtrl = $controller('BrowserCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
